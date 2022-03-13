#ifndef IPU_BASE_HPP
#define IPU_BASE_HPP

#include <poplar/Graph.hpp>
#include <poplar/Engine.hpp>

#include "swatlib/similarity.h"
#include "swatlib/encoding.h"
#include "ipu_config.h"

using namespace poplar;
namespace ipu {

void addCycleCount(Graph& graph, program::Sequence& mainProgram, const std::string handle);
uint64_t getTotalCycles(Engine& engine, const std::string handle);
double calculateGCUPS(uint64_t cellCount, double elapsedSeconds);

Type formatToType(const std::string& format);

TypeTraits typeToTrait(const Type& t);

void convertSimilarityMatrix(Target& target, const Type& t, swatlib::Matrix<int8_t> matrix, void** buffer);

inline void addCodelets(Graph& graph);

inline int extractScoreSW(Engine& engine, const std::string& sA, const std::string& sB);

class IPUAlgorithm {
    std::vector<poplar::Device> devices;
protected:

    int thread_id;
    int ipus;
    // std::unique_ptr<Engine> engine;
    // Engine* engine = nullptr;
    std::vector<poplar::Engine*> engines;
public:
    SWConfig config;
    IPUAlgorithm(SWConfig config, int thread_id, int ipu_count);

    Graph createGraph();

    // Needs to be called in child class
    void createEngine(Graph& graph, std::vector<program::Program> programs, std::string hash = "", bool use_cache = false);

    std::vector<poplar::Device>& getDevices();
    poplar::Graph getGraph();
    poplar::Target getTarget();

    double getTileClockFrequency();

    ~IPUAlgorithm();
};

}

#endif // IPU_BASE_HPP