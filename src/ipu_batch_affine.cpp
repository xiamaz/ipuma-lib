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

class FillCallback final : public poplar::StreamCallback {
 public:
  using Result = poplar::StreamCallback::Result;

  FillCallback(msd::channel<SubmittedBatch*>& value, std::map<slotToken, SubmittedBatch*>& results, std::mutex& rmutex, size_t size, size_t ipu_id) : ch(value), resultTable(results), size(size), id(ipu_id), result_mutex(rmutex) {}

  Result prefetch(void* __restrict p) noexcept override {
    PLOGE.printf("Enter PreFETCH, id %d", id);
    // We do this to handle closed streams;
    for (auto b : ch) {
      PLOGE.printf("Do PreFETCH, id %d", id);
      pushBatch(p, b);
      return Result::Success;
    }
    PLOGE.printf("Exit PreFETCH, id %d", id);
    if (ch.closed()) {
      close(p);
    }
    return Result::Success;
  }

  void fetch(void* __restrict p) noexcept override {
    PLOGE.printf("FETCH, %d", id);
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
    PLOGD.printf("Send teardown message to IPU, id %d", id);
    std::vector<int32_t> aa(size + 1, 0);
    memcpy(p, aa.data(), aa.size() * 4);
  }
  inline void pushBatch(void* __restrict p, SubmittedBatch* b) {
    // Wireformat: inputbuffer+slotToken
    // TODO!!!!!!!!!: Are the offsets correct?????
    size_t inputsize = abs((char*)b->inputs_end - (char*)b->inputs_begin);
    memcpy(p, (char*)b->inputs_begin, inputsize);
    int* ip = reinterpret_cast<int*>(&reinterpret_cast<char*>(p)[inputsize]);
    ip[0] = b->slot + 1;
    b->runTick.tick();
    result_mutex.lock();
    resultTable.insert({b->slot, b});
    result_mutex.unlock();
  }
  size_t size;
  size_t id;
  msd::channel<SubmittedBatch*>& ch;
  // TODO: Add mutex!
  // TODO: move this out of the hot worker?
  std::mutex &result_mutex;
  std::map<slotToken, SubmittedBatch*>& resultTable;
};

class RecvCallback final : public poplar::StreamCallback {
 public:
  using Result = poplar::StreamCallback::Result;

  RecvCallback(std::map<slotToken, SubmittedBatch*>& results, std::mutex& rmutex) : resultTable(results), result_mutex(rmutex) {}

  Result prefetch(void* __restrict p) noexcept override {
    PLOGE.printf("NOOOOOOOOOOOOOOOOOOOOOOOOOO");
    exit(1);
    return Result::NotAvailable;
  }

  void fetch(void* __restrict p) noexcept override {
    pullBatch(p);
  }

  void complete() noexcept override {}

 private:
  inline void pullBatch(void* __restrict p) {
    // Wireformat: slotToken+outbuffer
    int32_t* ip = reinterpret_cast<int32_t*>(p);
    slotToken st = ip[0] - 1;
    if (st == -1) {
      return;
    }
    result_mutex.lock();
    auto b = resultTable[st];
    *b->cyclesH2D = readTime(&ip[1]);
    *b->cyclesInner = readTime(&ip[1 + 2]);
    b->runTick.tock();
    // PLOGE.printf("Cycles: %lu, %lu", *b.cyclesH2D, *b.cyclesInner);
    resultTable.erase(st);
    result_mutex.unlock();
    memcpy((char*)b->results_begin, &ip[1 + 2 + 2], abs((char*)b->results_end - (char*)b->results_begin));
    b->signal_done->close();
  }

  uint64_t readTime(int32_t* mem) {
    uint32_t cycles[2];
    memcpy(cycles, mem, 4 * 2);
    uint64_t totalCycles = (((uint64_t)cycles[1]) << 32) | cycles[0];
    return totalCycles;
  }

  std::mutex &result_mutex;
  std::map<slotToken, SubmittedBatch*>& resultTable;
};

/**
 * Streamable IPU graph for SW
 */
std::vector<program::Program> buildGraph(Graph& graph, VertexType vtype, unsigned long activeTiles, unsigned long maxAB,
                                         unsigned long bufSize, unsigned long maxBatches,
                                         const swatlib::Matrix<int8_t> similarityData, int gapInit, int gapExt,
                                         bool use_remote_buffer, int transmissionPrograms, int buf_rows, bool forward_only, int ioTiles) {
  int buffers = transmissionPrograms;

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
    }

    TypeTraits traits = typeToTrait(sType);
    void* similarityBuffer;
    convertSimilarityMatrix(target, sType, similarityData, &similarityBuffer);
    Tensor similarity;
    similarity = graph.addConstant(sType, {m, n}, similarityBuffer, traits, false, "similarity");
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

    OptionFlags streamOptions({});
    program::Sequence undef;
    auto frontCs = graph.addComputeSet("SmithWaterman");
    for (int i = 0; i < activeTiles; ++i) {
      int tileIndex = i % tileCount;
      VertexRef vtx = graph.addVertex(frontCs, vertexTypeToIpuLabel(vtype),
                                      {
                                          {"bufSize", bufSize},
                                          {"maxAB", maxAB},
                                          {"gapInit", gapInit},
                                          {"gapExt", gapExt},
                                          {"maxNPerTile", maxBatches},
                                          {"Seqs", Seqs[i]},
                                          {"Meta", CompMeta[i]},
                                          {"simMatrix", similarity},
                                          {"score", Scores[i]},
                                          {"ARange", ARanges[i]},
                                          {"BRange", BRanges[i]},
                                          {"forwardOnly", forward_only},
                                      });
      // graph.setFieldSize(vtx["C"], maxAB * workerMultiplier);
      auto C_T = graph.addVariable(sType, {maxAB * workerMultiplier}, "C[" + std::to_string(i) + "]");
      graph.setTileMapping(C_T, tileIndex);
      graph.connect(vtx["C"], C_T);
      // undef.add(program::WriteUndef(C_T));

      // graph.setFieldSize(vtx["bG"], maxAB * workerMultiplier);
      auto bG_T = graph.addVariable(sType, {maxAB * workerMultiplier}, "bG[" + std::to_string(i) + "]");
      graph.setTileMapping(bG_T, tileIndex);
      graph.connect(vtx["bG"], bG_T);
      // undef.add(program::WriteUndef(bG_T));

      // if (vtype == VertexType::stripedasm) {
      //   assert(false && "Not Implemented");
      //   graph.connect(vtx["simWidth"], m);
      //   graph.setFieldSize(vtx["tS"], maxAB * workerMultiplier);
      // }
      // if (vtype == VertexType::multistriped || vtype == VertexType::multistripedasm) {
      //   assert(false && "Not Implemented");
      //   graph.setFieldSize(vtx["tS"], maxAB * target.getNumWorkerContexts());
      //   graph.setFieldSize(vtx["locks"], target.getNumWorkerContexts());
      // }
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
    if (use_remote_buffer) {
      const char* units[] = {"B", "kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
      int b_size = inputs_tensor.numElements() * 4;
      PLOGI.printf("Use a RemoteBuffer of bytes %.2f", (double)b_size);
      int i = 0;
      for (; b_size > 1024; i++) {
        b_size /= 1024;
      }
      PLOGI.printf("Use a RemoteBuffer of %dx(Banksize %.2f %s)", buf_rows, (double)b_size, units[i]);
      auto offset = graph.addVariable(INT, {1}, "Remote Buffer Offset");
      graph.setTileMapping(offset, 0);
      // graph.setInitialValue<int>(offset, {0});
      // auto cs = graph.addComputeSet("AddMod");
      // auto vtx = graph.addVertex(cs, "AddMod");
      // graph.connect(vtx["val"], offset[0]);
      // graph.setInitialValue<int>(vtx["mod"], buf_rows);
      // graph.setTileMapping(vtx, 0);

      // 256 MB row limit?
      // https://docs.graphcore.ai/projects/poplar-user-guide/en/latest/poplar_programs.html#remote-buffer-restrictions
      auto host_stream_offset = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT_N(buffer_share), INT, 1);
      h2d_prog_concat.add(poplar::program::Copy(host_stream_offset, offset));
      auto memid = REMOTE_MEMORY_N(buffer_share);
      auto remote_mem = graph.addRemoteBuffer(memid, inputs_tensor.elementType(), inputs_tensor.numElements(), buf_rows, true, true);
      // h2d_prog_concat.add(poplar::program::Copy(remote_mem, inputs_tensor, offset, "Copy Inputs from Remote->IPU"));
      // d2h_prog_concat.add(program::Execute(cs));

      device_stream_concat = graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL_N(buffer_share), INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements());
      d2h_prog_concat.add(poplar::program::Copy(outputs_tensor, device_stream_concat, "Copy Outputs from IPU->Host"));
    } else {
      auto host_stream_concat = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT_N(buffer_share), INT, inputs_tensor.numElements() + 1, ReplicatedStreamMode::REPLICATE, {{"splitLimit", std::to_string(264 * 1024 * 1024)}});
      device_stream_concat = graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL_N(buffer_share), INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements() + 1 + 2 + 2);
      auto inT = concat({inputs_tensor.flatten(), slotT.flatten()});
      PLOGE.printf("Input Buffer size = %lu bytes", inT.numElements() * 4);
      h2d_prog_concat.add(poplar::program::Copy(host_stream_concat, inT));
    }

    auto print_tensors_prog = program::Sequence({
        // program::PrintTensor("Alens", Alens),
        // program::PrintTensor("Blens", Blens),
        program::PrintTensor("CompMeta", CompMeta),
        // program::PrintTensor("Scores", Scores),
    });
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
  // return {prog, d2h_prog_concat, print_tensors_prog};
}

// TODO: Better default for work_queue!!!
SWAlgorithm::SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, size_t slotCap, size_t ipuCount)
    : IPUAlgorithm(config, thread_id, ipuCount), algoconfig(algoconfig), work_queue({10'000}), resultTable({}) {
  const auto totalComparisonsCount = algoconfig.getTotalNumberOfComparisons();
  slot_size = slotCap;
  slot_avail.resize(algoconfig.transmissionPrograms);
  std::fill(slot_avail.begin(), slot_avail.end(), slotCap);
  slots.resize(algoconfig.transmissionPrograms);
  for (size_t i = 0; i < algoconfig.transmissionPrograms; i++) {
    slots[i].resize(slot_size);
    std::fill(slots[i].begin(), slots[i].end(), false);
  }

  scores.resize(totalComparisonsCount);
  a_range_result.resize(totalComparisonsCount);
  b_range_result.resize(totalComparisonsCount);

  Graph graph = createGraph();

  auto similarityMatrix = swatlib::selectMatrix(config.similarity, config.matchValue, config.mismatchValue, config.ambiguityValue);
  std::vector<program::Program> programs =
      buildGraph(graph, algoconfig.vtype, algoconfig.tilesUsed, algoconfig.maxAB, algoconfig.bufsize, algoconfig.maxBatches,
                 similarityMatrix, config.gapInit, config.gapExtend, algoconfig.useRemoteBuffer, algoconfig.transmissionPrograms, slot_size, algoconfig.forwardOnly, algoconfig.ioTiles);

  std::hash<std::string> hasher;
  auto s = json{algoconfig, config};
  std::stringstream ss("");
  graph.outputComputeGraph(ss, programs);
  size_t hash = hasher(s.dump() + ss.str());
  createEngine(graph, programs, std::to_string(hash));
  run_executor();
}

SWAlgorithm::SWAlgorithm(ipu::SWConfig config, IPUAlgoConfig algoconfig)
    : SWAlgorithm::SWAlgorithm(config, algoconfig, 0) {}

BlockAlignmentResults SWAlgorithm::get_result() { return {scores, a_range_result, b_range_result}; }

void SWAlgorithm::checkSequenceSizes(const IPUAlgoConfig& algoconfig, const std::vector<int>& SeqSizes) {
  if (SeqSizes.size() > algoconfig.maxBatches * algoconfig.tilesUsed) {
    PLOGW << "Sequence has more elements than the maxBatchsize";
    PLOGW << "Seq.size() = " << SeqSizes.size();
    PLOGW << "max comparisons = " << algoconfig.tilesUsed << " * " << algoconfig.maxBatches;
    throw std::runtime_error("Input sequence (A) is over the max length.");
  }

  for (const auto& s : SeqSizes) {
    if (s > algoconfig.maxAB) {
      PLOGW << "Sequence size in seq " << s << " > " << algoconfig.maxAB;
      exit(1);
    }
  }
}

void SWAlgorithm::checkSequenceSizes(const IPUAlgoConfig& algoconfig, const RawSequences& Seqs) {
  std::vector<int> seqSizes;
  std::transform(Seqs.begin(), Seqs.end(), std::back_inserter(seqSizes), [](const auto& s) { return s.size(); });
  checkSequenceSizes(algoconfig, seqSizes);
}

void SWAlgorithm::fillMNBuckets(Algorithm algo, partition::BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int offset) {
  switch (algo) {
    case Algorithm::fillFirst:
      partition::fillFirst(map, Seqs, Cmps, offset);
      break;
    case Algorithm::roundRobin:
      partition::roundRobin(map, Seqs, Cmps, offset);
      break;
    case Algorithm::greedy:
      partition::greedy(map, Seqs, Cmps, offset);
      break;
  }
}

void SWAlgorithm::fillBuckets(Algorithm algo, partition::BucketMap& map, const RawSequences& A, const RawSequences& B, int offset) {
  switch (algo) {
    case Algorithm::fillFirst:
      partition::fillFirst(map, A, B, offset);
      break;
    case Algorithm::roundRobin:
      partition::roundRobin(map, A, B, offset);
      break;
    case Algorithm::greedy:
      partition::greedy(map, A, B, offset);
      break;
  }
}

size_t SWAlgorithm::getSeqsOffset(const IPUAlgoConfig& config) { return 0; }

size_t SWAlgorithm::getMetaOffset(const IPUAlgoConfig& config) {
  return config.getTotalBufsize32b();
}

void SWAlgorithm::refetch() {
  assert(false && "Not implemented");
  // engine->run(1);
}

std::tuple<int, slotToken> SWAlgorithm::unpack_slot(slotToken s) {
  return {s / slot_size, s % slot_size};
}

void SWAlgorithm::upload(int32_t* inputs_begin, int32_t* inputs_end, slotToken slot) {
  auto [lot, sid] = unpack_slot(slot);
  PLOGD.printf("Slot is %d, lot is %d", sid, lot);
  swatlib::TickTock rbt;
  rbt.tick();
  assert(false && "Not implemented");
  // engine->copyToRemoteBuffer(inputs_begin, REMOTE_MEMORY_N(lot), sid);
  rbt.tock();
  auto transferTime = rbt.duration<std::chrono::milliseconds>();
  size_t totalTransferSize = algoconfig.getInputBufferSize32b() * 4;
  auto transferBandwidth = (double)totalTransferSize / ((double)transferTime / 1000.0) / 1e6;
  PLOGI.printf("Transfer rate to remote buffer %.3f mb/s. %ld bytes in %.2fms", transferBandwidth, totalTransferSize, transferTime);
}

void Job::join() {
  for (auto done : *done_signal) {
    throw std::runtime_error("done_signal should be closed");
  }
  tick.tock();
}

Job* SWAlgorithm::async_submit_prepared_remote_compare(int32_t* inputs_begin, int32_t* inputs_end, int32_t* results_begin, int32_t* results_end) {
  auto job = new Job();
  job->done_signal = new msd::channel<int>();
  slotToken sid = (rand() % 100000000) + 1;
  uint64_t h2d_cycles, inner_cycles;
  job->sb = SubmittedBatch{inputs_begin, inputs_end, results_begin, results_end, sid, job->done_signal, &(job->h2dCycles), &(job->innerCycles), {}};
  job->tick.tick();
  &(job->sb) >> this->work_queue;
  return job;
}

void SWAlgorithm::blocking_join_prepared_remote_compare(Job& job) {
  job.join();
  // release_slot(slot_token);
  PLOGD << "Total engine run time (in s): " << static_cast<double>(job.tick.duration<std::chrono::milliseconds>()) / 1000.0;

#ifdef IPUMA_DEBUG
  // auto [lot, sid] = unpack_slot(slot_token);
  auto cyclesOuter = job.h2dCycles + job.innerCycles;
  auto cyclesInner = job.innerCycles;
  auto timeOuter = static_cast<double>(cyclesOuter) / getTarget().getTileClockFrequency();
  auto timeInner = static_cast<double>(cyclesInner) / getTarget().getTileClockFrequency();
  PLOGD << "Poplar cycle count: " << cyclesInner << "/" << cyclesOuter << " computed time (in s): " << timeInner << "/"
        << timeOuter;

  int buffers = algoconfig.transmissionPrograms;
  int ibufsize = (algoconfig.bufsize / 4) + (algoconfig.bufsize % 4 == 0 ? 0 : 1);
  auto scale_bufsize = [=](int bufsize, int lot) {
    return static_cast<size_t>((bufsize / buffers + (bufsize % buffers == 0 ? 0 : 1)) * (lot + 1));
  };
  // int32_t* meta_input = inputs_begin + scale_bufsize(getMetaOffset(algoconfig), lot);
  int32_t* meta_input = job.sb.inputs_begin + getMetaOffset(algoconfig);

  // GCUPS computation
  uint64_t cellCount = 0;
  uint64_t dataCount = 0;
  for (size_t i = 0; i < algoconfig.getTotalNumberOfComparisons(); i++) {
    auto a_len = meta_input[4 * i];
    auto b_len = meta_input[4 * i + 2];
    cellCount += a_len * b_len;
    dataCount += a_len + b_len;
    // PLOGW << a_len << " : blen : " << b_len;
  }

  double GCUPSOuter = static_cast<double>(cellCount) / timeOuter / 1e9;
  double GCUPSInner = static_cast<double>(cellCount) / timeInner / 1e9;
  PLOGD << "Poplar estimated cells(" << cellCount << ") GCUPS " << GCUPSInner << "/" << GCUPSOuter;

  // dataCount - actual data content transferred
  // totalTransferSize - size of buffer being transferred
  double totalTransferSize = algoconfig.getInputBufferSize32b() * 4;

  auto transferTime = timeOuter - timeInner;
  auto transferInfoRatio = static_cast<double>(dataCount) / totalTransferSize * 100;
  auto transferBandwidth = totalTransferSize / transferTime / 1e6;
  auto transferBandwidthPerVertex = transferBandwidth / algoconfig.tilesUsed;
  PLOGD << "Transfer time: " << transferTime << "s estimated bandwidth: " << transferBandwidth
        << "mb/s, per vertex: " << transferBandwidthPerVertex << "mb/s";
#endif
}

void SWAlgorithm::prepared_remote_compare(int32_t* inputs_begin, int32_t* inputs_end, int32_t* results_begin, int32_t* results_end, slotToken slot_token) {
  release_slot(slot_token);
  auto job = async_submit_prepared_remote_compare(inputs_begin, inputs_end, results_begin, results_end);
  blocking_join_prepared_remote_compare(*job);
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
  if (algoconfig.useRemoteBuffer) {
    slot = queue_slot(maxbucket);
    upload(&*inputs.begin() + 2, &*inputs.end() - 2, slot);
  }
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

std::string SWAlgorithm::printTensors() {
  // std::stringstream graphstream;
  // engine->setPrintTensorStream(graphstream);
  // engine->run(2);
  // return graphstream.str();
  assert(false && "Not implemented");
  return "";
}

void SWAlgorithm::transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end) {
  int numComparisons = mapping_end - mapping_begin;
  transferResults(results_begin, results_end, mapping_begin, mapping_end, scores_begin, scores_end, arange_begin, arange_end, brange_begin, brange_end, numComparisons);
}
void SWAlgorithm::transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end, int numComparisons) {
  size_t results_size = results_end - results_begin;
  size_t result_part_size = results_size / 3;
  if (results_size % 3 != 0) {
    std::cout << results_size << "\n";
    throw std::runtime_error("Results buffer not 3 aligned.");
  }
  auto* results_score = results_begin;
  auto* results_arange = results_begin + result_part_size;
  auto* results_brange = results_begin + result_part_size * 2;

  for (int i = 0; i < numComparisons; ++i) {
    auto mapped_i = mapping_begin[i];
    scores_begin[i] = results_score[mapped_i];
    arange_begin[i] = results_arange[mapped_i];
    brange_begin[i] = results_brange[mapped_i];
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
  if (algoconfig.useRemoteBuffer) {
    slot = queue_slot(maxBucket);
    upload(&*inputs.begin(), &*inputs.end(), slot);
  }
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

  const auto encodeTable = swatlib::getEncoder(dtype).getCodeTable();
  const size_t seqs_offset = getSeqsOffset(algoconfig);

  int8_t* seqs = (int8_t*)inputs_begin + seqs_offset;
  int32_t* meta = inputs_begin + meta_offset;

  // PLOGW << "======";
  // omp_set_num_threads(16);
  // PLOGW << omp_get_num_threads();
  // PLOGW << "======";
  // #pragma omp parallel for num_threads(32)
  // for (const auto& bucketMapping : map.buckets) {
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
#pragma omp simd
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

slotToken SWAlgorithm::queue_slot(int max_bucket_size) {
  assert(buf_has_capacity(max_bucket_size));
  auto si = calculate_slot_region_index(max_bucket_size);
  int s = -1;
  for (size_t i = 0; i < slots.size(); i++) {
    if (!slots[si][i]) {
      s = i;
      break;
    }
  }
  assert(s != -1);
  slots[si][s] = true;
  slot_avail[si]--;
  last_slot = s;
  return s + (si * slot_size);
}

void SWAlgorithm::release_slot(slotToken i) {
  auto [lot, sid] = unpack_slot(i);
  assert(slots[lot][sid] == true);
  slots[lot][sid] = false;
  slot_avail[lot]++;
}

bool SWAlgorithm::slot_available(int max_buffer_size) {
  return slot_avail[calculate_slot_region_index(max_buffer_size)] > 0;
}

int SWAlgorithm::calculate_slot_region_index(const IPUAlgoConfig& algoconfig, int max_buffer_size) {
  auto ret = (int)ceil((max_buffer_size / (double)algoconfig.bufsize) * static_cast<double>(algoconfig.transmissionPrograms)) - 1;
  PLOGD.printf("Calculated slot region %d from maxbuf %d and bufsize %d", ret, max_buffer_size, algoconfig.bufsize);
  return ret;
}

int SWAlgorithm::calculate_slot_region_index(int max_buffer_size) {
  auto ret = (int)ceil((max_buffer_size / (double)algoconfig.bufsize) * static_cast<double>(algoconfig.transmissionPrograms)) - 1;
  PLOGD.printf("Calculated slot region %d from maxbuf %d and bufsize %d", ret, max_buffer_size, algoconfig.bufsize);
  return ret;
}

void SWAlgorithm::run_executor() {
  // Connect output
  for (size_t i = 0; i < ipus; i++) {
    std::unique_ptr<RecvCallback> rb{new RecvCallback(resultTable, tableMutex)};
    engines[i]->connectStreamToCallback(STREAM_CONCAT_ALL_N(0), std::move(rb));

    // Connect input
    std::unique_ptr<FillCallback> cb{new FillCallback(work_queue, resultTable, tableMutex,algoconfig.getInputBufferSize32b(), i)};
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