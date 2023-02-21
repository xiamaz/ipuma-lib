#include "ipu_batch_affine.h"

#include <omp.h>
#include <plog/Log.h>

#include <cmath>
#include <iostream>
#include <nlohmann/json.hpp>
#include <poplar/CycleCount.hpp>
#include <poplar/Graph.hpp>
#include <poplar/SyncType.hpp>
#include <popops/Zero.hpp>
#include <poputil/TileMapping.hpp>
#include <string>
#include <thread>

#include "msd/channel.hpp"
#include "swatlib/swatlib.h"

namespace ipu {
namespace batchaffine {

using json = nlohmann::json;

#define STREAM_CONCAT_ALL_N(n) "concat-read-all" + std::to_string(n)
#define HOST_STREAM_CONCAT_N(n) "host-stream-concat" + std::to_string(n)
#define REMOTE_MEMORY_N(n) "inputs-remotebuffer" + std::to_string(n)

#define CYCLE_COUNT_OUTER_N(n) "cycle-count-outer" + std::to_string(n)
#define CYCLE_COUNT_INNER_N(n) "cycle-count-inner" + std::to_string(n)

/**
 * @brief Callback for pushing jobs to the IPU.
 * 
 */
class FillCallback final : public poplar::StreamCallback {
 public:
  using Result = poplar::StreamCallback::Result;

  FillCallback(
    msd::channel<Job*>& value,
    JobMap& results,
    std::mutex& rmutex,
    size_t size,
    size_t ipu_id
  ) : ch(value), resultTable(results), size(size), id(ipu_id), result_mutex(rmutex) {}

  Result prefetch(void* __restrict p) noexcept override {
    PLOGV.printf("Enter PreFETCH, ipuId %d", id);
    // We do this to handle closed streams;
    for (auto b : ch) {
      PLOGV.printf("Do PreFETCH, ipuId %d", id);
      pushBatch(p, b);
      return Result::Success;
    }
    PLOGV.printf("Exit PreFETCH, ipuId %d", id);
    if (ch.closed()) {
      close(p);
    }
    return Result::NotAvailable;
  }

  void fetch(void* __restrict p) noexcept override {
    PLOGV.printf("FETCH, ipuId %d", id);
    // We do this to handle closed streams;
    for (auto b : ch) {
      pushBatch(p, b);
      return;
    }
    if (ch.closed()) {
      close(p);
    }
    return;
  }

  void complete() noexcept override {}

 private:
  void close(void* __restrict p) {
    PLOGW.printf("Send teardown message to IPU, ipuId %d", id);
    std::vector<int32_t> aa(size + 1, 0);
    memcpy(p, aa.data(), aa.size() * 4);
  }
  void pushBatch(void* __restrict p, Job* j) {
    // Wireformat: inputbuffer+JobId
    const auto& inputBuffer = j->batch->inputs;
    assert(inputBuffer.size() > 0);
    size_t inputSize = inputBuffer.size() * sizeof(inputBuffer[0]);

    PLOGD << "Input buffer size " << inputSize << " defined size is " << size;
    memcpy(p, inputBuffer.data(), inputSize);

    // insert job ID at the end of the transfer buffer
    int* ip = reinterpret_cast<int*>(&reinterpret_cast<char*>(p)[inputSize]);
    ip[0] = j->id;

    j->batch->tick.tick();
    PLOGV << "PushBatch slot " << j->id;

    // inserting job into resultTable so that we can add result later
    result_mutex.lock();
    resultTable.insert({j->id, j});
    result_mutex.unlock();
  }
  size_t size;
  size_t id;
  msd::channel<Job*>& ch;
  // TODO: Add mutex!
  // TODO: move this out of the hot worker?
  std::mutex &result_mutex;
  JobMap& resultTable;
};

/**
 * @brief Callback for getting back results from the IPU
 * 
 */
class RecvCallback final : public poplar::StreamCallback {
 public:
  using Result = poplar::StreamCallback::Result;

  RecvCallback(JobMap& results, std::mutex& rmutex) : resultTable(results), result_mutex(rmutex) {}

  Result prefetch(void* __restrict p) noexcept override {
    PLOGF.printf("NOOOOOOOOOOOOOOOOOOOOOOOOOO");
    exit(1);
    return Result::NotAvailable;
  }

  void fetch(void* __restrict p) noexcept override {
    pullBatch(p);
  }

  void complete() noexcept override {}

 private:
  void pullBatch(void* __restrict p) {
    // Wireformat: slotToken+outbuffer
    int32_t* ip = reinterpret_cast<int32_t*>(p);
    JobId jobId = ip[0];
    PLOGV << "PullBatch JobId " << jobId;
    if (jobId == -1 || jobId == 0) {
      return;
    }

    // getting Job object based on id
    result_mutex.lock();
    auto* j = resultTable[jobId];
    j->h2dCycles = readTime(&ip[1]);
    j->innerCycles = readTime(&ip[1 + 2]);
    j->batch->tick.tock();
    resultTable.erase(jobId);
    result_mutex.unlock();

    auto& resultBuffer = j->batch->results;
    size_t bufferSize = resultBuffer.size() * sizeof(resultBuffer[0]);
    memcpy((char*)resultBuffer.data(), &ip[1 + 2 + 2], bufferSize);
    j->done_signal->close();

    PLOGV << "PullBatch EXIT jobId " << jobId;
  }

  uint64_t readTime(int32_t* mem) {
    uint32_t cycles[2];
    memcpy(cycles, mem, 4 * 2);
    uint64_t totalCycles = (((uint64_t)cycles[1]) << 32) | cycles[0];
    return totalCycles;
  }

  std::mutex& result_mutex;
  JobMap& resultTable;
};

void Job::join() {
  for (auto done : *done_signal) {
    throw std::runtime_error("done_signal should be closed");
  }
  tick.tock();
}

/**
 * Streamable IPU graph for SW
 */
std::vector<program::Program> buildGraph(Graph& graph, VertexType vtype, unsigned long activeTiles, unsigned long maxSequenceLength,
                                         unsigned long bufSize, unsigned long maxComparisonsPerVertex,
                                         const swatlib::Matrix<int8_t> similarityData, int gapInit, int gapExt,
                                         bool forward_only, int ioTiles,
                                         int xDrop, double bandPercentageXDrop, int seedLength, int metaSize
                                         ) {
  // TODO: remove transmissionPrograms
  int buffers = 1;

  std::vector<program::Program> progs(buffers);
  for (size_t buffer_share = 0; buffer_share < buffers; buffer_share++) {
    program::Sequence prog;
    auto target = graph.getTarget();
    int tileCount = target.getTilesPerIPU();

    int ibufsize = (bufSize / 4) + (bufSize % 4 == 0 ? 0 : 1);
    size_t seq_scaled_size = static_cast<size_t>((ibufsize / buffers + (ibufsize % buffers == 0 ? 0 : 1)) * (buffer_share + 1));
    PLOGD.printf("Set seq size %d of total %d", seq_scaled_size, bufSize / 4);
    Tensor Seqs = graph.addVariable(INT, {activeTiles, seq_scaled_size}, "Seqs");

    // Metadata structure
    Tensor CompMeta = graph.addVariable(INT, {activeTiles, maxComparisonsPerVertex * (metaSize)}, "CompMeta");

    auto [m, n] = similarityData.shape();

    poplar::Type sType;
    int workerMultiplier = 1;
    switch (vtype) {
      case VertexType::cpp:
        sType = INT;
        break;
      case VertexType::assembly:
        sType = FLOAT;
        break;
      case VertexType::multi:
        sType = INT;
        workerMultiplier = target.getNumWorkerContexts();
        break;
      case VertexType::multiasm:
        sType = FLOAT;
        workerMultiplier = target.getNumWorkerContexts();
        break;
      case VertexType::xdrop:
        sType = INT; 
        break;
      case VertexType::multixdrop:
        // This ok?
        sType = INT; 
        workerMultiplier = target.getNumWorkerContexts();
        break;
      case VertexType::multibandxdrop:
        // This ok?
        sType = INT; 
        workerMultiplier = target.getNumWorkerContexts();
        break;
      case VertexType::xdropseedextend:
        // This ok?
        sType = INT; 
        workerMultiplier = target.getNumWorkerContexts();
        break;
      case VertexType::greedyxdrop:
        // This ok?
        sType = INT; 
        break;
      default:
        PLOGF.printf("Unknown vtype %d", vtype);
        throw "Unknown vtype";
    }

    TypeTraits traits = typeToTrait(sType);
    void* similarityBuffer;
    convertSimilarityMatrix(target, sType, similarityData, &similarityBuffer);
    Tensor similarity = graph.addConstant(sType, {m, n}, similarityBuffer, traits, false, "similarity");
    free(similarityBuffer);

    // TypeTraits traits = typeToTrait(sType);
    // void* similarityBuffer;
    // convertSimilarityMatrix(target, sType, similarityData, &similarityBuffer);
    // Tensor similarity;
    // if (vtype == VertexType::stripedasm) {
    //   similarity = graph.addConstant(sType, {m * n}, similarityBuffer, traits, false, "similarity");
    // } else {
    //   similarity = graph.addConstant(sType, {m, n}, similarityBuffer, traits, false, "similarity");
    // }
    // free(similarityBuffer);

    Tensor Scores = graph.addVariable(INT, {activeTiles, maxComparisonsPerVertex}, "Scores");
    Tensor ARanges = graph.addVariable(INT, {activeTiles, maxComparisonsPerVertex}, "ARanges");
    Tensor BRanges = graph.addVariable(INT, {activeTiles, maxComparisonsPerVertex}, "BRanges");

    graph.setTileMapping(similarity, 0);
    if (ioTiles != 0) {
      PLOGW.printf("Use %d I/O Tiles", ioTiles);
      for (int i = 0; i < activeTiles; ++i) {
        int tileIndex = activeTiles + std::max((int)ceil(ioTiles * (double)i / (double)activeTiles) - 1, 0);
        // PLOGD.printf("I/O Tile %d", tileIndex);
        graph.setTileMapping(Seqs[i], tileIndex);
        graph.setTileMapping(CompMeta[i], tileIndex);
      }

      for (int i = 0; i < activeTiles; ++i) {
        int tileIndex = i % tileCount;
        graph.setTileMapping(Scores[i], tileIndex);
        graph.setTileMapping(ARanges[i], tileIndex);
        graph.setTileMapping(BRanges[i], tileIndex);
      }
    } else {
      for (int i = 0; i < activeTiles; ++i) {
        int tileIndex = i % tileCount;
        graph.setTileMapping(Seqs[i], tileIndex);
        graph.setTileMapping(CompMeta[i], tileIndex);

        graph.setTileMapping(Scores[i], tileIndex);
        graph.setTileMapping(ARanges[i], tileIndex);
        graph.setTileMapping(BRanges[i], tileIndex);
      }
    }

    auto label = vertexTypeToIpuLabel(vtype);
    if (vtype == VertexType::multibandxdrop) {
      label += "<" + std::to_string(xDrop) + ">";
    } else if (vtype == VertexType::multixdrop) {
      label += "<" + std::to_string(xDrop) + ">"; 
    } else if (vtype == VertexType::xdropseedextend) {
      label += "<" + std::to_string(xDrop) + ">";
    }
    PLOGD.printf("Use Vertex class: %s", label.c_str());
    OptionFlags streamOptions({});
    program::Sequence undef;
    auto frontCs = graph.addComputeSet("SmithWaterman");
    for (int i = 0; i < activeTiles; ++i) {
      int tileIndex = i % tileCount;
      VertexRef vtx = graph.addVertex(frontCs, label,
                                      {
                                          {"bufSize", bufSize},
                                          {"gapInit", gapInit},
                                          {"gapExt", gapExt},
                                          {"maxNPerTile", maxComparisonsPerVertex},
                                          {"Seqs", Seqs[i]},
                                          {"Meta", CompMeta[i]},
                                          {"score", Scores[i]},
                                      });

      if (vtype == VertexType::xdrop) {
        auto k_T = graph.addVariable(sType, {3, (maxSequenceLength+2) * workerMultiplier}, "K[" + std::to_string(i) + "]");
        graph.setTileMapping(k_T, tileIndex);
        graph.connect(vtx["maxSequenceLength"], maxSequenceLength);
        graph.connect(vtx["K1"], k_T[0]);
        graph.connect(vtx["K2"], k_T[1]);
        graph.connect(vtx["K3"], k_T[2]);
        graph.connect(vtx["simMatrix"], similarity);
        PLOGE << "OUTDATED";
      } else if (vtype == VertexType::multixdrop) {
        auto k_T = graph.addVariable(sType, {2, (maxSequenceLength+2) * workerMultiplier}, "K[" + std::to_string(i) + "]");
        graph.connect(vtx["maxSequenceLength"], maxSequenceLength);
        graph.setTileMapping(k_T, tileIndex);
        graph.connect(vtx["K1"], k_T[0]);
        graph.connect(vtx["K2"], k_T[1]);
        graph.connect(vtx["simMatrix"], similarity);
      } else if (vtype == VertexType::xdropseedextend) {
        auto k_T = graph.addVariable(sType, {2, (maxSequenceLength+2) * workerMultiplier}, "K[" + std::to_string(i) + "]");
        graph.connect(vtx["maxSequenceLength"], maxSequenceLength);
        graph.setTileMapping(k_T, tileIndex);
        graph.connect(vtx["K1"], k_T[0]);
        graph.connect(vtx["K2"], k_T[1]);
        graph.connect(vtx["simMatrix"], similarity);
        graph.connect(vtx["seedLength"], seedLength);
      } else if (vtype == VertexType::multibandxdrop) {
        int scaledMaxAB = maxSequenceLength * bandPercentageXDrop;
        auto k_T = graph.addVariable(sType, {2, ((size_t)scaledMaxAB+2+2) * (size_t) workerMultiplier}, "K[" + std::to_string(i) + "]");
        graph.connect(vtx["maxSequenceLength"], scaledMaxAB);
        graph.setTileMapping(k_T, tileIndex);
        graph.connect(vtx["K1"], k_T[0]);
        graph.connect(vtx["K2"], k_T[1]);
        graph.connect(vtx["simMatrix"], similarity);
        PLOGE << "OUTDATED";
      } else if (vtype == VertexType::greedyxdrop) {
        int mis = -2;
        int mat = 2;
        int X = 10;
        int xdrop_offset = ((X + mat / 2) / (mat - mis)) + 1;

        graph.connect(vtx["maxSequenceLength"], maxSequenceLength);
        Tensor R0 = graph.addVariable(INT, {2 * maxSequenceLength});
        graph.setTileMapping(R0, tileIndex);
        Tensor R1 = graph.addVariable(INT, {2 * maxSequenceLength});
        graph.setTileMapping(R1, tileIndex);
        Tensor T = graph.addVariable(INT, {maxSequenceLength+maxSequenceLength+xdrop_offset+1});
        graph.setTileMapping(T, tileIndex);
        graph.connect(vtx["VTR0"], R0);
        graph.connect(vtx["VTR1"], R1);
        graph.connect(vtx["VTT"], T);
      } else {
        graph.connect(vtx["maxSequenceLength"], maxSequenceLength);
        graph.connect(vtx["ARange"], ARanges[i]);
        graph.connect(vtx["BRange"], BRanges[i]);
        graph.connect(vtx["forwardOnly"], forward_only);
        // graph.setFieldSize(vtx["C"], maxSequenceLength * workerMultiplier);
        auto C_T = graph.addVariable(sType, {maxSequenceLength * workerMultiplier}, "C[" + std::to_string(i) + "]");
        graph.setTileMapping(C_T, tileIndex);
        graph.connect(vtx["C"], C_T);
        // undef.add(program::WriteUndef(C_T));

        // graph.setFieldSize(vtx["bG"], maxSequenceLength * workerMultiplier);
        auto bG_T = graph.addVariable(sType, {maxSequenceLength * workerMultiplier}, "bG[" + std::to_string(i) + "]");
        graph.setTileMapping(bG_T, tileIndex);
        graph.connect(vtx["bG"], bG_T);
        graph.connect(vtx["simMatrix"], similarity);
        // undef.add(program::WriteUndef(bG_T));
      }
      graph.setTileMapping(vtx, tileIndex);
      graph.setPerfEstimate(vtx, 1);
    }

    PLOGD.printf("Seqs offset %d in bytes %d", Seqs.flatten().numElements(), Seqs.flatten().numElements() * 4);
    auto inputs_tensor = concat({Seqs.flatten(), CompMeta.flatten()});
    auto outputs_tensor = concat({Scores.flatten(), ARanges.flatten(), BRanges.flatten()});

    program::Sequence h2d_prog_concat;
    program::Sequence d2h_prog_concat;

    auto slotT = graph.addVariable(INT, {1}, "SlotToken+1");
    graph.setTileMapping(slotT, tileCount - 1);
    DataStream device_stream_concat;

    auto host_stream_concat = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT_N(buffer_share), INT, inputs_tensor.numElements() + 1, ReplicatedStreamMode::REPLICATE, {{"splitLimit", std::to_string(264 * 1024 * 1024)}, {"bufferingDepth", "1"}});
    device_stream_concat = graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL_N(buffer_share), INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements() + 1 + 2 + 2);
    auto inT = concat({inputs_tensor.flatten(), slotT.flatten()});
    PLOGE.printf("Input Buffer size = %lu bytes", inT.numElements() * 4);
    h2d_prog_concat.add(poplar::program::Copy(host_stream_concat, inT));

    program::Sequence main_prog;
    main_prog.add(program::Execute(frontCs));
    // #ifdef IPUMA_DEBUG
    //    addCycleCount(graph, main_prog, CYCLE_COUNT_INNER_N(buffer_share));
    // #endif
    Tensor cyclesInner = poplar::cycleCount(graph, main_prog, 0, SyncType::EXTERNAL);
    Tensor cyclesH2D = poplar::cycleCount(graph, h2d_prog_concat, 0, SyncType::EXTERNAL);
    prog.add(h2d_prog_concat);
    prog.add(main_prog);
    // prog.add(program::PrintTensor(cycles));
    auto outT = concat({slotT.flatten(), cyclesH2D.reinterpret(INT).flatten(), cyclesInner.reinterpret(INT).flatten(), outputs_tensor.flatten()});
    PLOGE.printf("Output Buffer size = %lu bytes", outT.numElements() * 4);
    d2h_prog_concat.add(poplar::program::Copy(outT, device_stream_concat));
    prog.add(d2h_prog_concat);
    // #ifdef IPUMA_DEBUG
    //     addCycleCount(graph, prog, CYCLE_COUNT_OUTER_N(buffer_share));
    // #endif
    program::Sequence reps;
    program::Sequence test;
    reps.add(prog);
    reps.add(program::RepeatWhileTrue(test, slotT[0], prog));
    // reps.add( program::PrintTensor("EXIT AS SLOT IS",slotT));
    // reps.add(program::Repeat(10000, prog));
    progs[buffer_share] = reps;
  }
  return progs;
}

// TODO: Better default for work_queue!!!
SWAlgorithm::SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, size_t ipuCount, bool runExecutor)
    : IPUAlgorithm(config, thread_id, ipuCount), algoconfig(algoconfig), work_queue({WORK_QUEUE_SIZE}), resultTable({}) {
  Graph graph = createGraph();

  auto similarityMatrix = swatlib::selectMatrix(config.similarity, config.matchValue, config.mismatchValue, config.ambiguityValue);
  std::vector<program::Program> programs = buildGraph(
    graph,
    algoconfig.vtype,
    algoconfig.numVertices,
    algoconfig.maxSequenceLength,
    algoconfig.vertexBufferSize,
    algoconfig.maxComparisonsPerVertex,
    similarityMatrix,
    config.gapInit,
    config.gapExtend,
    algoconfig.forwardOnly,
    algoconfig.ioTiles,
    algoconfig.xDrop,
    algoconfig.bandPercentageXDrop,
    algoconfig.seedLength,
    algoconfig.getMetaStructSize32b());

  std::hash<std::string> hasher;
  auto s = json{algoconfig, config};
  std::stringstream ss("");
  graph.outputComputeGraph(ss, programs);
  size_t hash = hasher(s.dump() + ss.str());
  createEngine(graph, programs, std::to_string(hash));

  if (runExecutor) run_executor();
}

SWAlgorithm::SWAlgorithm(ipu::SWConfig config, IPUAlgoConfig algoconfig)
    : SWAlgorithm::SWAlgorithm(config, algoconfig, 0) {}

void SWAlgorithm::checkSequenceSizes(const IPUAlgoConfig& algoconfig, const RawSequences& Seqs) {
  if (Seqs.size() > algoconfig.maxComparisonsPerVertex * algoconfig.numVertices) {
    PLOGW << "Sequence has more elements than the maxBatchsize";
    PLOGW << "Seq.size() = " << Seqs.size();
    PLOGW << "max comparisons = " << algoconfig.numVertices << " * " << algoconfig.maxComparisonsPerVertex;
    throw std::runtime_error("Input sequence (A) is over the max length.");
  }

  std::vector<int> seqSizes;
  for (const auto& sequence : Seqs) {
    const auto size = sequence.size();
    if (size > algoconfig.maxSequenceLength) {
      PLOGW << "Sequence size in seq " << size << " > " << algoconfig.maxSequenceLength;
      exit(1);
    }
  }
}

Job* SWAlgorithm::async_submit(Batch* batch) {
  if (batch != nullptr) {
    auto job = new Job();
    job->done_signal = new msd::channel<int>();
    job->id = (rand() % 100000000) + 1;
    job->batch = batch;
    job->tick.tick();
    job >> this->work_queue;
    return job;
  } else {
    PLOGE << "null pointer passed as batch.";
  }
  return nullptr;
}

void computeJobMetrics(const Job& job, double tileFrequency, size_t numberVertices, size_t maxComparisonsPerVertex) {
  auto cyclesOuter = job.h2dCycles + job.innerCycles;
  auto cyclesInner = job.innerCycles;
  auto timeOuter = static_cast<double>(cyclesOuter) / tileFrequency;
  auto timeInner = static_cast<double>(cyclesInner) / tileFrequency;

  // GCUPS computation
  double GCUPSOuter = static_cast<double>(job.batch->cellCount) / timeOuter / 1e9;
  double GCUPSInner = static_cast<double>(job.batch->cellCount) / timeInner / 1e9;

  // transfer computation
  double totalTransferSize = job.batch->inputs.size() * sizeof(job.batch->inputs[0]);
  auto transferTime = static_cast<double>(job.h2dCycles) / tileFrequency;
  auto transferInfoRatio = static_cast<double>(job.batch->dataCount) / totalTransferSize * 100;
  auto transferBandwidth = totalTransferSize / transferTime / 1e6;
  auto transferBandwidthPerVertex = transferBandwidth / numberVertices;

  double timeJob = job.tick.seconds();
  double timeBatch = job.batch->tick.seconds();

  json logData = {
    {"job_id", job.id},
    {"cycles_h2d", job.h2dCycles},
    {"cycles_inner", job.innerCycles},
    {"cycles_outer", cyclesOuter},
    {"time_job", timeJob},
    {"time_batch", timeBatch},
    {"time_inner", timeInner},
    {"time_outer", timeOuter},
    {"time_transfer", transferTime },
    {"gcups_inner", GCUPSInner},
    {"gcups_outer", GCUPSOuter},
    {"transfer_bandwidth", transferBandwidth},
    {"transfer_bandwidth_per_vertex", transferBandwidthPerVertex},
    {"transfer_info_ratio", transferInfoRatio},
    {"cell_count", job.batch->cellCount},
    {"comparison_occupancy", job.batch->numComparisons / (double) (numberVertices * maxComparisonsPerVertex) * 100.0},
  };

  PLOGD << "JOBLOG: " << logData.dump();
}

void SWAlgorithm::blocking_join(Job& job) {
  job.join();
  // release_slot(slot_token);

#ifdef IPUMA_DEBUG
  computeJobMetrics(job, getTileClockFrequency(), algoconfig.numVertices, algoconfig.maxComparisonsPerVertex);
#endif
}

std::vector<Batch> SWAlgorithm::create_batches(const RawSequences& seqs, const Comparisons& cmps) {
  auto batches = ipu::create_batches(seqs, cmps, algoconfig, config);
  return batches;
}

void SWAlgorithm::run_executor() {
  // Connect output
  for (size_t i = 0; i < ipus; i++) {
    std::unique_ptr<RecvCallback> rb{new RecvCallback(resultTable, tableMutex)};
    engines[i]->connectStreamToCallback(STREAM_CONCAT_ALL_N(0), std::move(rb));

    // Connect input
    std::unique_ptr<FillCallback> cb{new FillCallback(work_queue, resultTable, tableMutex, algoconfig.getInputBufferSize32b(), i)};
    engines[i]->connectStreamToCallback(HOST_STREAM_CONCAT_N(0), std::move(cb));

    // Server
    executor_procs.push_back(std::thread([](poplar::Engine* engine, int ipu_id) {
      PLOGI.printf("Run Engine");
      if (engine == nullptr) {
        PLOGE.printf("Run Engine is null");
      }
      engine->run(0);
      PLOGW.printf("Engine exited with IPU_ID %d", ipu_id);
    }, engines[i], i));
  }
}

SWAlgorithm::~SWAlgorithm() {
  PLOGD.printf("Initiaite close");
  work_queue.close();
  PLOGD.printf("Initiate done");
  for (auto&& ep : executor_procs) {
    ep.join();
  }
  PLOGD.printf("Join done");
}

}  // namespace batchaffine
}  // namespace ipu