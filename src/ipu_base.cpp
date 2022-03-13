#include <plog/Log.h>

#include <linux/limits.h>
#include <string>
#include <stdexcept>

#include <libgen.h>
#include <unistd.h>

#include <poplar/Engine.hpp>
#include <poplar/SyncType.hpp>
#include <poplar/CycleCount.hpp>
#include <poplar/Program.hpp>
#include <poplar/DeviceManager.hpp>
#include <poplar/IPUModel.hpp>
#include <poplar/VariableMappingMethod.hpp>

#include <poputil/TileMapping.hpp>

#include "popops/ElementWise.hpp"
#include "popops/codelets.hpp"

#include "swatlib/swatlib.h"

#include "ipu_base.h"

using namespace poplar;
namespace ipu {

void addCycleCount(Graph& graph, program::Sequence& mainProgram, const std::string handle) {
    Tensor cycles = cycleCount(graph, mainProgram, 1, SyncType::EXTERNAL);
    graph.createHostRead(handle, cycles);
}

uint64_t getTotalCycles(Engine& engine, const std::string handle) {
  uint32_t cycles[2];
  engine.readTensor(handle, &cycles, &cycles + 1);
  uint64_t totalCycles = (((uint64_t)cycles[1]) << 32) | cycles[0];
  return totalCycles;
}

double calculateGCUPS(uint64_t cellCount, double elapsedSeconds) {
    return static_cast<double>(cellCount) / elapsedSeconds;
}

TypeTraits typeToTrait(const Type& t) {
    if (t == INT) {
        return TypeTraits::make<int32_t>();
    } else if (t == FLOAT) {
        return TypeTraits::make<float>();
    } else if (t == HALF) {
        return TypeTraits::make<IeeeHalf>();
    } else if (t == SHORT) {
        return TypeTraits::make<int16_t>();
    } else {
        throw std::runtime_error("Unsupported type for trait conversion: " + t.toString());
    }
}

void convertSimilarityMatrix(Target& target, const Type& t, swatlib::Matrix<int8_t> matrix, void** buffer) {
    auto [m, n] = matrix.shape();
    void* b;
    int8_t* matrixp = matrix.data();
    if (t == INT) {
        b = malloc(m * n * sizeof(int32_t));
        auto s = static_cast<int32_t*>(b);
        for (int i = 0; i < m * n; ++i) {
            s[i] = static_cast<int32_t>(matrixp[i]);
        }
    } else if (t == FLOAT) {
        b = malloc(m * n * sizeof(float));
        auto s = static_cast<float*>(b);
        for (int i = 0; i < m * n; ++i) {
            s[i] = static_cast<float>(matrixp[i]);
        }
    } else if (t == HALF) {
        b = malloc(m * n * 2);
        float* fbuf = (float*) malloc(m * n * sizeof(float));
        for (int i = 0; i < m * n; ++i) {
            fbuf[i] = static_cast<float>(matrixp[i]);
        }
        copyFloatToDeviceHalf(target, fbuf, b, m * n);
        free(fbuf);
    } else if (t == SHORT) {
        b = malloc(m * n * sizeof(int16_t));
        auto s = static_cast<int16_t*>(b);
        for (int i = 0; i < m * n; ++i) {
            s[i] = static_cast<int16_t>(matrixp[i]);
        }
    } else {
        throw std::runtime_error("Unknown type format: " + t.toString());
    }
    *buffer = b;
}

std::string get_selfpath() {
    char buff[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
      buff[len] = '\0';
      return std::string(buff);
    }
    return "";
}

inline void addCodelets(Graph& graph) {
    auto selfPath = get_selfpath();
    auto rootPath = std::string(dirname(dirname(&selfPath[0])));
    PLOGD << "Loading codelets from "<< rootPath << "/bin/codelets/algoipu.gp";
    graph.addCodelets(rootPath + "/bin/codelets/algoipu.gp");
}

inline int extractScoreSW(Engine& engine, const std::string& sA, const std::string& sB) {
    unsigned long m = sA.length() + 1;
    unsigned long n = sB.length() + 1;

    swatlib::Matrix<short> S(m, n);
    swatlib::Matrix<char> D(m, n);

    engine.readTensor("s-read", &*S.begin(), &*S.end());
    engine.readTensor("d-read", &*D.begin(), &*D.end());

    auto [x, y] = S.argmax();
    return S(x, y);
}

IPUAlgorithm::IPUAlgorithm(SWConfig config, int thread_id, int ipu_count) : config(config), thread_id(thread_id), ipus(ipu_count) {
    devices.resize(ipus);
    engines.resize(ipus, nullptr);
    auto manager = poplar::DeviceManager::createDeviceManager();
    // Attempt to attach to a single IPU:

    auto devs = manager.getDevices(TargetType::IPU, 1);
    auto attached = 0;
    for (size_t i = 0; i < devs.size(); i++) {
        auto& d = devs.at(thread_id + i);
        PLOGD.printf("Trying to aquired IPU ID %d", i); 
        if (d.attach()) {
            devices[attached] = std::move(d);
            attached++;
            PLOGI.printf("Aquired IPU ID %d", i); 
        }         
        if (attached >= ipus) break;
    }
    if (attached < ipus) {
        throw std::runtime_error("Could not aquire enough devices.");
    }
}

Graph IPUAlgorithm::createGraph() {
    auto target = getTarget();
    Graph graph(target);
    addCodelets(graph);
    popops::addCodelets(graph);
    return graph;
}

// Needs to be called in child class
void IPUAlgorithm::createEngine(Graph& graph, std::vector<program::Program> programs, const std::string hash, bool use_cache) {
    auto& device = getDevices();
    poplar::OptionFlags engineOptions;
    engineOptions.set("exchange.streamBufferOverlap", "none");
    engineOptions.set("exchange.enablePrefetch", "true");
    engineOptions.set("streamCallbacks.multiThreadMode", "collaborative");
    engineOptions.set("streamCallbacks.numWorkerThreads", "2");
    engineOptions.set("target.hostSyncTimeout", "0");
    // engineOptions.set("exchange.streamBufferOverlap", "none");
    // engineOptions.set("exchange.enablePrefetch", "true");
    std::fstream infile; 
    auto filename = "./"+hash+".poplar_exec";
    infile.open(filename);
    Executable executable;
    if (infile.good() && use_cache) {
        PLOGW.printf("Load Executeable %s.", filename.c_str());
        executable = Executable::deserialize(infile);
        infile.close();
    } else {
        executable = poplar::compileGraph(graph, programs, engineOptions);
        std::ofstream execFileCbor;
        PLOGW.printf("Write Executeable %s.", filename.c_str());
        execFileCbor.open(filename);
        executable.serialize(execFileCbor);
        execFileCbor.close();
    }

    for (size_t i = 0; i < ipus; i++) {
        std::fstream xinfile; 
        xinfile.open(filename);
        PLOGW.printf("Load Executeable %s.", filename.c_str());
        executable = Executable::deserialize(xinfile);
        xinfile.close();
        engines[i] = new Engine(std::move(executable), engineOptions);
        engines[i]->load(devices[i]);
    }
}

double IPUAlgorithm::getTileClockFrequency() {
    return getTarget().getTileClockFrequency();
}

IPUAlgorithm::~IPUAlgorithm() {
    for (size_t i = 0; i < ipus; i++) {
        delete engines[i];
        devices[i].detach();
    }
}

poplar::Target IPUAlgorithm::getTarget() {
    return devices[0].getTarget();
}

std::vector<poplar::Device>& IPUAlgorithm::getDevices() {
    return devices;
}

poplar::Graph IPUAlgorithm::getGraph() {
    return std::move(poplar::Graph(getTarget()));
}

}