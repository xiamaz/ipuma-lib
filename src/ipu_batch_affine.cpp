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
    PLOGV.printf("Enter PreFETCH, id %d", id);
    // We do this to handle closed streams;
    for (auto b : ch) {
      PLOGV.printf("Do PreFETCH, id %d", id);
      pushBatch(p, b);
      return Result::Success;
    }
    PLOGV.printf("Exit PreFETCH, id %d", id);
    if (ch.closed()) {
      close(p);
    }
    return Result::NotAvailable;
  }

  void fetch(void* __restrict p) noexcept override {
    PLOGV.printf("FETCH, %d", id);
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
    PLOGW.printf("Send teardown message to IPU, id %d", id);
    std::vector<int32_t> aa(size + 1, 0);
    memcpy(p, aa.data(), aa.size() * 4);
  }
  void pushBatch(void* __restrict p, Job* j) {
    // Wireformat: inputbuffer+JobId
    const auto& inputBuffer = j->batch->inputs;
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
    if (jobId == -1) {
      return;
    }
    PLOGV << "PullBatch JobId " << jobId;

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

BlockAlignmentResults Batch::get_result() {
  std::vector<int32_t> scores(numComparisons);
  std::vector<int32_t> a_range_result(numComparisons);
  std::vector<int32_t> b_range_result(numComparisons);

  // TODO this should not be hardcoded here but passed as a config to make layout changes easier.
  int aOffset = results.size() / 3;
  int bOffset = aOffset * 2;

  // mapping is now responsibility of the caller
  for (int i = 0; i < numComparisons; ++i) {
    int aindex = i + aOffset;
    int bindex = i + bOffset;
    scores[i] = results[i];
    a_range_result[i] = results[aindex];
    b_range_result[i] = results[bindex];
  }

  return {scores, a_range_result, b_range_result};
}

void Job::join() {
  for (auto done : *done_signal) {
    throw std::runtime_error("done_signal should be closed");
  }
  tick.tock();
}

/**
 * Streamable IPU graph for SW
 */
std::vector<program::Program> buildGraph(Graph& graph, VertexType vtype, unsigned long activeTiles, unsigned long maxAB,
                                         unsigned long bufSize, unsigned long maxBatches,
                                         const swatlib::Matrix<int8_t> similarityData, int gapInit, int gapExt,
                                         bool forward_only, int ioTiles,
                                         int xDrop, double bandPercentageXDrop
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
    Tensor CompMeta = graph.addVariable(INT, {activeTiles, maxBatches * 4}, "CompMeta");

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
      case VertexType::greedyxdrop:
        // This ok?
        sType = INT; 
        break;
      default:
        PLOGF.printf("Unknown vtype $d", vtype);
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

    Tensor Scores = graph.addVariable(INT, {activeTiles, maxBatches}, "Scores");
    Tensor ARanges = graph.addVariable(INT, {activeTiles, maxBatches}, "ARanges");
    Tensor BRanges = graph.addVariable(INT, {activeTiles, maxBatches}, "BRanges");

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
                                          {"maxNPerTile", maxBatches},
                                          {"Seqs", Seqs[i]},
                                          {"Meta", CompMeta[i]},
                                          {"score", Scores[i]},
                                      });

      if (vtype == VertexType::xdrop) {
        auto k_T = graph.addVariable(sType, {3, (maxAB+2) * workerMultiplier}, "K[" + std::to_string(i) + "]");
        graph.setTileMapping(k_T, tileIndex);
        graph.connect(vtx["maxAB"], maxAB);
        graph.connect(vtx["K1"], k_T[0]);
        graph.connect(vtx["K2"], k_T[1]);
        graph.connect(vtx["K3"], k_T[2]);
        graph.connect(vtx["simMatrix"], similarity);
      } else if (vtype == VertexType::multixdrop) {
        auto k_T = graph.addVariable(sType, {2, (maxAB+2) * workerMultiplier}, "K[" + std::to_string(i) + "]");
        graph.connect(vtx["maxAB"], maxAB);
        graph.setTileMapping(k_T, tileIndex);
        graph.connect(vtx["K1"], k_T[0]);
        graph.connect(vtx["K2"], k_T[1]);
        graph.connect(vtx["simMatrix"], similarity);
      } else if (vtype == VertexType::multibandxdrop) {
        int scaledMaxAB = maxAB * bandPercentageXDrop;
        auto k_T = graph.addVariable(sType, {2, ((size_t)scaledMaxAB+2+2) * (size_t) workerMultiplier}, "K[" + std::to_string(i) + "]");
        graph.connect(vtx["maxAB"], scaledMaxAB);
        graph.setTileMapping(k_T, tileIndex);
        graph.connect(vtx["K1"], k_T[0]);
        graph.connect(vtx["K2"], k_T[1]);
        graph.connect(vtx["simMatrix"], similarity);
      } else if (vtype == VertexType::greedyxdrop) {
        int mis = -2;
        int mat = 2;
        int X = 10;
        int xdrop_offset = ((X + mat / 2) / (mat - mis)) + 1;

        graph.connect(vtx["maxAB"], maxAB);
        Tensor R0 = graph.addVariable(INT, {2 * maxAB});
        graph.setTileMapping(R0, tileIndex);
        Tensor R1 = graph.addVariable(INT, {2 * maxAB});
        graph.setTileMapping(R1, tileIndex);
        Tensor T = graph.addVariable(INT, {maxAB+maxAB+xdrop_offset+1});
        graph.setTileMapping(T, tileIndex);
        graph.connect(vtx["VTR0"], R0);
        graph.connect(vtx["VTR1"], R1);
        graph.connect(vtx["VTT"], T);
      } else {
        graph.connect(vtx["maxAB"], maxAB);
        graph.connect(vtx["ARange"], ARanges[i]);
        graph.connect(vtx["BRange"], BRanges[i]);
        graph.connect(vtx["forwardOnly"], forward_only);
        // graph.setFieldSize(vtx["C"], maxAB * workerMultiplier);
        auto C_T = graph.addVariable(sType, {maxAB * workerMultiplier}, "C[" + std::to_string(i) + "]");
        graph.setTileMapping(C_T, tileIndex);
        graph.connect(vtx["C"], C_T);
        // undef.add(program::WriteUndef(C_T));

        // graph.setFieldSize(vtx["bG"], maxAB * workerMultiplier);
        auto bG_T = graph.addVariable(sType, {maxAB * workerMultiplier}, "bG[" + std::to_string(i) + "]");
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
    algoconfig.tilesUsed,
    algoconfig.maxAB,
    algoconfig.bufsize,
    algoconfig.maxBatches,
    similarityMatrix,
    config.gapInit,
    config.gapExtend,
    algoconfig.forwardOnly,
    algoconfig.ioTiles,
    algoconfig.xDrop,
    algoconfig.bandPercentageXDrop
  );

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
  if (Seqs.size() > algoconfig.maxBatches * algoconfig.tilesUsed) {
    PLOGW << "Sequence has more elements than the maxBatchsize";
    PLOGW << "Seq.size() = " << Seqs.size();
    PLOGW << "max comparisons = " << algoconfig.tilesUsed << " * " << algoconfig.maxBatches;
    throw std::runtime_error("Input sequence (A) is over the max length.");
  }

  std::vector<int> seqSizes;
  for (const auto& sequence : Seqs) {
    const auto size = sequence.size();
    if (size > algoconfig.maxAB) {
      PLOGW << "Sequence size in seq " << size << " > " << algoconfig.maxAB;
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

void computeJobMetrics(const Job& job, double tileFrequency, size_t numberVertices) {
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

  double timeJob = static_cast<double>(job.tick.duration<std::chrono::milliseconds>()) / 1000.0;
  double timeBatch = static_cast<double>(job.batch->tick.duration<std::chrono::milliseconds>()) / 1000.0;

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
  };

  PLOGD << "JOBLOG: " << logData.dump();
}

void SWAlgorithm::blocking_join(Job& job) {
  job.join();
  // release_slot(slot_token);

#ifdef IPUMA_DEBUG
  computeJobMetrics(job, getTileClockFrequency(), algoconfig.tilesUsed);
#endif
}

std::vector<Batch> SWAlgorithm::create_batches(const RawSequences& seqs, const Comparisons& cmps) {
  partition::BucketMap map(
    algoconfig.tilesUsed,
    algoconfig.maxBatches,
    algoconfig.bufsize
  );
  fillBuckets(algoconfig.fillAlgo, map, seqs, cmps, 0);
  std::vector<Batch> batches;
  const auto inputBufferSize = algoconfig.getInputBufferSize32b();
  const auto mappingBufferSize = algoconfig.getTotalNumberOfComparisons();

  batches.push_back({});

  // fill data into batches

  json logData = {
    {"batches_created", batches.size()}
  };
  PLOGD << "BATCHCREATE: " << logData;
  return batches;
}

void SWAlgorithm::compare_local_many(const std::vector<std::string>& A, const std::vector<std::string>& B) {
  swatlib::TickTock tPrepare, tCompare;

  std::vector<int32_t*> mapping;
  std::vector<int32_t*> inputs;

  tPrepare.tick();
  prepare_local_many(config, algoconfig, A, B, inputs, mapping);
  tPrepare.tock();
  PLOGE << "MAPPING SIZE " << mapping.size();
  PLOGE << "MAPPING FIRST VAL " << ((int32_t*)mapping[0])[0];

  PLOGE << "Starting comparisons";

  for (int batch = 0; batch < mapping.size(); ++batch) {
    int32_t* input_begin = inputs[batch];
    int32_t* input_end = input_begin + algoconfig.getInputBufferSize32b();
    size_t results_size = scores.size() + a_range_result.size() + b_range_result.size();
    std::vector<int32_t> results(results_size);
    tCompare.tick();
    PLOGE << "prepared_remote_compare";
    prepared_remote_compare(input_begin, input_end, &*results.begin(), &*results.end());
    tCompare.tock();
    int32_t* mapping_begin = mapping[batch];
    int32_t* mapping_end = mapping_begin + algoconfig.getTotalNumberOfComparisons();
    PLOGE << "transfer results";
    transferResults(&*results.begin(), &*results.end(), mapping_begin, mapping_end, &*scores.begin(), &*scores.end(),
                    &*a_range_result.begin(), &*a_range_result.end(), &*b_range_result.begin(), &*b_range_result.end());
  }
  PLOGE << "Comps Done";

  double prepare_time = static_cast<double>(tPrepare.duration<std::chrono::milliseconds>()) / 1e3;
  double compare_time = static_cast<double>(tCompare.duration<std::chrono::milliseconds>()) / 1e3;
  double total_time = prepare_time + compare_time;
  double prepare_perc = prepare_time / total_time * 100;
  PLOGD << "Total time: " << total_time << " prepare(" << prepare_time << ") compare(" << compare_time << ") preprocessing is "
        << prepare_perc << "% of total";

  PLOGE << "Exit func";
}

void SWAlgorithm::compare_local(const std::vector<std::string>& A, const std::vector<std::string>& B, bool errcheck) {
  swatlib::TickTock tPrepare, tCompare;
  std::vector<int> mapping(A.size(), 0);
  size_t inputs_size = algoconfig.getInputBufferSize32b();
  std::vector<int32_t> inputs(inputs_size + 4);

  size_t results_size = scores.size() + a_range_result.size() + b_range_result.size();
  std::vector<int32_t> results(results_size + 4);

  inputs[0] = 0xDEADBEEF;
  inputs[1] = 0xDEADBEEF;
  inputs[inputs_size + (1) + 1] = 0xFEEBDAED;
  inputs[inputs_size + (1) + 2] = 0xFEEBDAED;
  tPrepare.tick();
  int maxbucket = prepare_remote(config, algoconfig, A, B, &*inputs.begin() + 2, &*inputs.end() - 2, mapping.data());
  tPrepare.tock();

  if (inputs[0] != 0xDEADBEEF || inputs[1] != 0xDEADBEEF) {
    std::vector<int32_t> subset(inputs.begin(), inputs.begin() + 10);
    PLOGW << "Canary begin overwritten " << swatlib::printVector(subset);
    exit(1);
  }
  if (inputs[inputs_size + 2] != 0xFEEBDAED || inputs[inputs_size + 3] != 0xFEEBDAED) {
    std::vector<int32_t> subset(inputs.end() - 10, inputs.end());
    PLOGW << "Canary end overwritten " << swatlib::printVector(subset);
    exit(1);
  }

  // std::vector<int32_t> unord_scores(scores), unord_mismatches(mismatches), unord_a_range(a_range_result),
  // unord_b_range(b_range_result);
  results[0] = 0xDEADBEEF;
  results[1] = 0xDEADBEEF;
  results[results_size + (1) + 1] = 0xFEEBDAED;
  results[results_size + (1) + 2] = 0xFEEBDAED;
  // prepared_remote_compare(a.data(), a_len.data(), b.data(), b_len.data(), unord_scores.data(), unord_mismatches.data(),
  // unord_a_range.data(),
  //                         unord_b_range.data());
  slotToken slot = 0;
  slot = queue_slot(maxbucket);
  tCompare.tick();
  prepared_remote_compare(&*inputs.begin() + 2, &*inputs.end() - 2, &*results.begin() + 2, &*results.end() - 2, slot);
  tCompare.tock();

  double prepare_time = static_cast<double>(tPrepare.duration<std::chrono::milliseconds>()) / 1e3;
  double compare_time = static_cast<double>(tCompare.duration<std::chrono::milliseconds>()) / 1e3;
  double total_time = prepare_time + compare_time;
  double prepare_perc = prepare_time / total_time * 100;
  PLOGD << "Total time: " << total_time << " prepare(" << prepare_time << ") compare(" << compare_time << ") preprocessing is "
        << prepare_perc << "% of total";

  // check canaries
  if (results[0] != 0xDEADBEEF || results[1] != 0xDEADBEEF) {
    std::vector<int32_t> subset(results.begin(), results.begin() + 10);
    PLOGW << "Canary begin overwritten " << swatlib::printVector(subset);
    exit(1);
  }
  if (results[results_size + 2] != 0xFEEBDAED || results[results_size + 3] != 0xFEEBDAED) {
    std::vector<int32_t> subset(results.end() - 10, results.end());
    PLOGW << "Canary end overwritten " << swatlib::printVector(subset);
    exit(1);
  }

  // PLOGD << "Unordered results: " << swatlib::printVector(std::vector<int32_t>(results.begin() + 2, results.begin() +
  // scores.size() + 2));

  // reorder results based on mapping
  int nthTry = 0;
  int sc;
retry:
  nthTry++;
  sc = 0;
  for (size_t i = 0; i < mapping.size(); ++i) {
    size_t mapped_i = mapping[i];
    scores[i] = results[mapped_i + 2];
    if (errcheck) {
      if (scores[i] >= algoconfig.maxAB) {
        PLOGW << "ERROR Expected " << A.size() << " valid comparisons. But got " << i << " instead.";
        PLOGW.printf("ERROR Thread %d received wrong data FIRST, try again data=%d, map_translate=%d\n", thread_id, scores[i],
                     mapping[i] + 2);
        refetch();
        goto retry;
      }
    }
    sc += scores[i] > 0;
    size_t a_range_offset = scores.size();
    size_t b_range_offset = a_range_offset + a_range_result.size();
    a_range_result[i] = results[a_range_offset + mapped_i + 2];
    b_range_result[i] = results[b_range_offset + mapped_i + 2];
  }
  if (errcheck) {
    if ((double)sc / A.size() < 0.5) {
      PLOGW << "ERROR Too many scores are 0, retry number " << (nthTry - 1);
      refetch();
      goto retry;
    }
  }
}

void SWAlgorithm::compare_mn_local(const std::vector<std::string>& Seqs, const Comparisons& Cmps, bool errcheck) {
  size_t inputs_size = algoconfig.getInputBufferSize32b();
  std::vector<int32_t> inputs(inputs_size);
  size_t results_size = scores.size() + a_range_result.size() + b_range_result.size();
  std::vector<int32_t> results(results_size);

  std::vector<int> mapping(Cmps.size(), 0);
  partition::BucketMap map(algoconfig.tilesUsed, algoconfig.maxBatches, algoconfig.bufsize);
  fillMNBuckets(algoconfig.fillAlgo, map, Seqs, Cmps, 0);
  fill_input_buffer(map, config.datatype, algoconfig, Seqs, Cmps, &*inputs.begin(), &*inputs.end(), mapping.data());
  // std::cout << "Mapping: " << swatlib::printVector(mapping) << "\n";
  int maxBucket = 0;
  for (const auto& bucket : map.buckets) {
    maxBucket = std::max(maxBucket, bucket.seqSize);
  }

#ifdef IPUMA_DEBUG
  int emptyBuckets = 0;
  long dataCount = 0;
  std::vector<int> bucketCmps;
  std::map<int, int> occurence;
  for (const auto& bucket : map.buckets) {
    if (bucket.cmps.size() == 0) emptyBuckets++;
    occurence[bucket.cmps.size()]++;
    bucketCmps.push_back(bucket.cmps.size());
    dataCount += bucket.seqSize;
  }
  std::stringstream ss;
  ss << "Map[";
  for (auto [k, v] : occurence) {
    ss << k << ": " << v << ",";
  }
  ss << "]";
  PLOGD << "Total number of buckets: " << map.numBuckets << " empty buckets: " << emptyBuckets;
  PLOGD << "Bucket size occurence: " << ss.str();
  int adjusted_bufsize = ((algoconfig.bufsize / algoconfig.transmissionPrograms) * (1 + calculate_slot_region_index(algoconfig, maxBucket)));
  double bucketPerc = static_cast<double>(maxBucket) / static_cast<double>(algoconfig.bufsize) * 100.0;
  double adjusted_bucketPerc = static_cast<double>(maxBucket) / static_cast<double>(adjusted_bufsize) * 100.0;
  PLOGD << "Max bucket: " << maxBucket << "/" << algoconfig.bufsize << " (" << bucketPerc << "%), adjusted " << maxBucket << "/" << adjusted_bufsize << " (" << adjusted_bucketPerc << "%)";
  double totalTransferSize = ((algoconfig.getInputBufferSize32b() * 4) / algoconfig.transmissionPrograms) * (1 + calculate_slot_region_index(algoconfig, maxBucket));
  auto transferInfoRatio = static_cast<double>(dataCount) / totalTransferSize * 100;
  PLOGD << "Transfer info/total: " << dataCount << "/" << totalTransferSize << " (" << transferInfoRatio << "%)";
#endif

  slotToken slot = 0;
  slot = queue_slot(maxBucket);
  prepared_remote_compare(&*inputs.begin(), &*inputs.end(), &*results.begin(), &*results.end(), slot);

  transferResults(&*results.begin(), &*results.end(), &*mapping.begin(), &*mapping.end(), &*scores.begin(), &*scores.end(),
                  &*a_range_result.begin(), &*a_range_result.end(), &*b_range_result.begin(), &*b_range_result.end());
}

void SWAlgorithm::fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& Seqs, const Comparisons& Cmps, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping) {
  auto encodeTable = swatlib::getEncoder(dtype).getCodeTable();
  const size_t seqs_offset = getSeqsOffset(algoconfig);
  size_t meta_offset = getMetaOffset(algoconfig);

  int maxBucket = 0;
  for (const auto& bucketMapping : map.buckets) {
    maxBucket = std::max(maxBucket, bucketMapping.seqSize);
  }

  int buffers = algoconfig.transmissionPrograms;
  int buffer_share = calculate_slot_region_index(algoconfig, maxBucket);
  int ibufsize = (algoconfig.bufsize / 4) + (algoconfig.bufsize % 4 == 0 ? 0 : 1);
  auto scale_bufsize = [=](int bufsize) {
    return static_cast<size_t>((bufsize / buffers + (bufsize % buffers == 0 ? 0 : 1)) * (buffer_share + 1));
  };
  size_t seq_scaled_size = scale_bufsize(ibufsize);

  PLOGD.printf("Detected max buffer of %d, this corresponds to lot %d, %d/%d", maxBucket, buffer_share, seq_scaled_size, ibufsize);
  PLOGD.printf("Old buffer %d/%d recalculated offsets %d/%d", seq_scaled_size, ibufsize, seq_scaled_size * algoconfig.tilesUsed, meta_offset);
  meta_offset = seq_scaled_size * algoconfig.tilesUsed;

  int8_t* seqs = (int8_t*)inputs_begin + seqs_offset;
  int32_t* meta = inputs_begin + meta_offset;

  for (int zzz = 0; zzz < map.buckets.size(); ++zzz) {
    const auto& bucketMapping = map.buckets[zzz];
    // const auto& bucketMapping : map.buckets) {
    const size_t offsetBuffer = bucketMapping.bucketIndex * scale_bufsize(algoconfig.getBufsize32b() * 4);
    const size_t offsetMeta = bucketMapping.bucketIndex * algoconfig.maxBatches * 4;
    auto* bucket_meta = meta + offsetMeta;
    auto* bucket_seq = seqs + offsetBuffer;

    for (const auto& sm : bucketMapping.seqs) {
      const char* seq;
      int seqSize;
      switch (sm.origin) {
        case partition::SequenceOrigin::A:
          throw std::runtime_error("Using A/B mapping with unordered comparison.");
          break;
        case partition::SequenceOrigin::B:
          throw std::runtime_error("Using A/B mapping with unordered comparison.");
          break;
        case partition::SequenceOrigin::unordered:
          seq = Seqs[sm.index].data();
          seqSize = Seqs[sm.index].size();
          break;
      }

      for (int j = 0; j < seqSize; ++j) {
        bucket_seq[sm.offset + j] = encodeTable[seq[j]];
      }
    }
    for (int i = 0; i < bucketMapping.cmps.size(); ++i) {
      const auto& cmpMapping = bucketMapping.cmps[i];
      bucket_meta[i * 4] = cmpMapping.sizeA;
      bucket_meta[i * 4 + 1] = cmpMapping.offsetA;
      bucket_meta[i * 4 + 2] = cmpMapping.sizeB;
      bucket_meta[i * 4 + 3] = cmpMapping.offsetB;

      mapping[cmpMapping.comparisonIndex] = map.cmpCapacity * bucketMapping.bucketIndex + i;
    }
  }
}

void SWAlgorithm::fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& A, const RawSequences& B, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping) {
  int maxBucket = 0;
  for (const auto& bucketMapping : map.buckets) {
    maxBucket = std::max(maxBucket, bucketMapping.seqSize);
  }

  int buffers = algoconfig.transmissionPrograms;
  int buffer_share = calculate_slot_region_index(algoconfig, maxBucket);
  int ibufsize = (algoconfig.bufsize / 4) + (algoconfig.bufsize % 4 == 0 ? 0 : 1);
  auto scale_bufsize = [=](int bufsize) {
    return static_cast<size_t>((bufsize / buffers + (bufsize % buffers == 0 ? 0 : 1)) * (buffer_share + 1));
  };
  size_t seq_scaled_size = scale_bufsize(ibufsize);

  PLOGD.printf("Detected max buffer of %d, this corresponds to lot %d, %d/%d", maxBucket, buffer_share, seq_scaled_size, ibufsize);
  size_t meta_offset = getMetaOffset(algoconfig);
  PLOGD.printf("Old buffer %d/%d recalculated offsets %d/%d", seq_scaled_size, ibufsize, seq_scaled_size * algoconfig.tilesUsed, meta_offset);
  meta_offset = seq_scaled_size * algoconfig.tilesUsed;

  size_t input_elems = inputs_end - inputs_begin;
  memset(inputs_begin, 0, input_elems * sizeof(int32_t));

  PLOGI << "Number of input elements " << input_elems;

  const auto encodeTable = swatlib::getEncoder(dtype).getCodeTable();
  const size_t seqs_offset = getSeqsOffset(algoconfig);

  int8_t* seqs = (int8_t*)inputs_begin + seqs_offset;
  int32_t* meta = inputs_begin + meta_offset;

 for (int zzz = 0; zzz < map.buckets.size(); ++zzz) {
    const auto& bucketMapping = map.buckets[zzz];
    const size_t offsetBuffer = bucketMapping.bucketIndex * scale_bufsize(algoconfig.getBufsize32b() * 4);
    const size_t offsetMeta = bucketMapping.bucketIndex * algoconfig.maxBatches * 4;
    auto* bucket_meta = meta + offsetMeta;
    auto* bucket_seq = seqs + offsetBuffer;

    for (const auto& sm : bucketMapping.seqs) {
      const char* seq;
      int seqSize;
      switch (sm.origin) {
        case partition::SequenceOrigin::A:
          seq = A[sm.index].data();
          seqSize = A[sm.index].size();
          break;
        case partition::SequenceOrigin::B:
          seq = B[sm.index].data();
          seqSize = B[sm.index].size();
          break;
        case partition::SequenceOrigin::unordered:
          throw std::runtime_error("Using unordered mapping with A/B comparison.");
          break;
      }
      // PLOGD << "Sequence size written " << seqSize;
#pragma omp simd
      for (int j = 0; j < seqSize; ++j) {
        bucket_seq[sm.offset + j] = encodeTable[seq[j]];
      }
      // PLOGD << "bucket_seq " << seqSize;
    }
    // PLOGD << "Here X1";
    for (int i = 0; i < bucketMapping.cmps.size(); ++i) {
      const auto& cmpMapping = bucketMapping.cmps[i];
      bucket_meta[i * 4] = cmpMapping.sizeA;
      bucket_meta[i * 4 + 1] = cmpMapping.offsetA;
      bucket_meta[i * 4 + 2] = cmpMapping.sizeB;
      bucket_meta[i * 4 + 3] = cmpMapping.offsetB;

      // PLOGW << "Here set mapping idx " << cmpMapping.comparisonIndex << " -> " <<  (map.cmpCapacity * bucketMapping.bucketIndex + i);
      // TODO: This is wrong!!!
      mapping[map.cmpCapacity * bucketMapping.bucketIndex + i] = cmpMapping.comparisonIndex;
      // This should be correct but is not.
      // mapping[cmpMapping.comparisonIndex] = map.cmpCapacity * bucketMapping.bucketIndex + i;
    }
    // PLOGD << "Here X2";
  }
    // PLOGD << "Here return";
}

void SWAlgorithm::prepare_local_many(
  const SWConfig& swconfig,
  const IPUAlgoConfig& algoconfig,
  const std::vector<std::string>& A,
  const std::vector<std::string>& B,
  std::vector<int32_t*>& inputs_begins,
  std::vector<int32_t*>& seqMappings
  ) {
  partition::BucketMap map(algoconfig.tilesUsed, algoconfig.maxBatches, algoconfig.bufsize);
  fillBuckets(algoconfig.fillAlgo, map, A, B, 0);


  size_t batches = map.numBuckets / algoconfig.tilesUsed;
  PLOGW << "NUM BATCHES!!!!!!!!!!!!!!!!!!!!!!!!: " << batches;

  // const auto resultBufferSize = algoconfig.getTotalNumberOfComparisons() * 3;

  inputs_begins.resize(batches);
  seqMappings.resize(batches);
  for (int i = 0; i<batches; i++) {
    inputs_begins[i] = (int32_t*) malloc(inputBufferSize * sizeof(int32_t));
    seqMappings[i] = (int32_t*) malloc(mappingBufferSize * sizeof(int32_t));
    PLOGW << "FRESH MAPPING VALUE " << seqMappings[0][0];

    partition::BucketMap maptmp(algoconfig.tilesUsed, algoconfig.maxBatches, algoconfig.bufsize);
    std::copy(map.buckets.begin() + i * algoconfig.tilesUsed, map.buckets.begin() + (i + 1 * algoconfig.tilesUsed), maptmp.buckets.begin());
    PLOGD << "Input buffer size " << inputBufferSize;

    fill_input_buffer(maptmp, swconfig.datatype, algoconfig, A, B, inputs_begins[i], inputs_begins[i]+inputBufferSize, seqMappings[i]);
  }
  PLOGD << "Finished prepare_local";
}

int SWAlgorithm::prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int32_t* inputs_begin, int32_t* inputs_end, int* seqMapping) {
  swatlib::TickTock preprocessTimer;
  std::vector<swatlib::TickTock> stageTimers(3);
  preprocessTimer.tick();
  checkSequenceSizes(algoconfig, A);
  // checkSequenceSizes(algoconfig, B);

  stageTimers[1].tick();
  partition::BucketMap map(algoconfig.tilesUsed, algoconfig.maxBatches, algoconfig.bufsize);
  fillBuckets(algoconfig.fillAlgo, map, A, B, 0);
  stageTimers[1].tock();

  stageTimers[2].tick();
  fill_input_buffer(map, swconfig.datatype, algoconfig, A, B, inputs_begin, inputs_end, seqMapping);
  stageTimers[2].tock();
  preprocessTimer.tock();

  PLOGD << "Total preprocessing time (in s): "
        << static_cast<double>(preprocessTimer.duration<std::chrono::milliseconds>()) / 1000.0;

  PLOGD << "Stage timers:";
  for (int i = 0; i < stageTimers.size(); ++i) {
    // PLOGD << i << ": " << stageTimers[i].duration<std::chrono::milliseconds>();
    PLOGD << i << ": " << stageTimers[i].accumulate_microseconds() / 1e3;
  }

  int maxBucket = 0;
  for (const auto& bucket : map.buckets) {
    maxBucket = std::max(maxBucket, bucket.seqSize);
  }
#ifdef IPUMA_DEBUG
  int emptyBuckets = 0;
  long dataCount = 0;
  std::vector<int> bucketCmps;
  std::map<int, int> occurence;
  for (const auto& bucket : map.buckets) {
    if (bucket.cmps.size() == 0) emptyBuckets++;
    occurence[bucket.cmps.size()]++;
    bucketCmps.push_back(bucket.cmps.size());
    dataCount += bucket.seqSize;
  }
  std::stringstream ss;
  ss << "Map[";
  for (auto [k, v] : occurence) {
    ss << k << ": " << v << ",";
  }
  ss << "]";
  PLOGD << "Total number of buckets: " << map.numBuckets << " empty buckets: " << emptyBuckets;
  PLOGD << "Bucket size occurence: " << ss.str();
  int adjusted_bufsize = ((algoconfig.bufsize / algoconfig.transmissionPrograms) * (1 + calculate_slot_region_index(algoconfig, maxBucket)));
  double bucketPerc = static_cast<double>(maxBucket) / static_cast<double>(algoconfig.bufsize) * 100.0;
  double adjusted_bucketPerc = static_cast<double>(maxBucket) / static_cast<double>(adjusted_bufsize) * 100.0;
  PLOGD << "Max bucket: " << maxBucket << "/" << algoconfig.bufsize << " (" << bucketPerc << "%), adjusted " << maxBucket << "/" << adjusted_bufsize << " (" << adjusted_bucketPerc << "%)";
  double totalTransferSize = ((algoconfig.getInputBufferSize32b() * 4) / algoconfig.transmissionPrograms) * (1 + calculate_slot_region_index(algoconfig, maxBucket));
  auto transferInfoRatio = static_cast<double>(dataCount) / totalTransferSize * 100;
  PLOGD << "Transfer info/total: " << dataCount << "/" << totalTransferSize << " (" << transferInfoRatio << "%)";
#endif

  return maxBucket;
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