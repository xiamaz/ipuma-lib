#include <string>
#include <cmath>
#include <iostream>

#include <plog/Log.h>

#include <poplar/Graph.hpp>
#include <poputil/TileMapping.hpp>
#include <popops/Zero.hpp>

#include "ipu_base.h"
#include "partition.h"
#include "swatlib/swatlib.h"
#include "ipu_batch_affine.h"

namespace ipu {
namespace batchaffine {

static const std::string STREAM_CONCAT_ALL = "concat-read-all";
static const std::string HOST_STREAM_CONCAT = "host-stream-concat";
static const std::string REMOTE_MEMORY = "inputs-remotebuffer";

static const std::string CYCLE_COUNT_OUTER = "cycle-count-outer";
static const std::string CYCLE_COUNT_INNER = "cycle-count-inner";

void printBucketMapping(BucketMapping bm) {
  std::cout << "Bucket:" << "\n";
  for (auto& c : bm.comparisons) {
    printf("%d: %d %d\n", c.cmpIndex, c.indexA, c.indexB);
  }
  std::cout << "=======================\n";
}

long long getCellCount(const std::vector<std::string>& A, const std::vector<std::string>& B) {
  long long cellCount = 0;
  if (A.size() != B.size()) {
    PLOGW << "Mismatch between size of A " << A.size() << " and size of B " << B.size();
  }
  // count cells based on 1:1 comparisons
  for (int i = 0; i < A.size(); ++i) {
    cellCount += A[i].size() * B[i].size();
  }
  return cellCount;
}

int IPUAlgoConfig::getBufsize32b() const {return std::ceil(static_cast<double>(bufsize) / 4.0); }

int IPUAlgoConfig::getTotalNumberOfComparisons() const { return tilesUsed * maxBatches; }

int IPUAlgoConfig::getMetaBufferSize32b() const {return getTotalNumberOfComparisons() * 4;};

int IPUAlgoConfig::getSequenceBufferSize8b() const { return tilesUsed * bufsize; }
int IPUAlgoConfig::getSequenceBufferSize32b() const { return tilesUsed * getBufsize32b(); }

int IPUAlgoConfig::getInputBufferSize32b() const { return getSequenceBufferSize32b() + getMetaBufferSize32b(); }

std::string vertexTypeToString(VertexType v) { return typeString[static_cast<int>(v)]; }

/**
 * Streamable IPU graph for SW
 */
std::vector<program::Program> buildGraph(Graph& graph, VertexType vtype, unsigned long activeTiles, unsigned long maxAB,
                                         unsigned long bufSize, unsigned long maxBatches,
                                         const swatlib::Matrix<int8_t> similarityData, int gapInit, int gapExt, bool use_remote_buffer, int buf_rows) {
  program::Sequence prog;

  auto target = graph.getTarget();
  int tileCount = target.getTilesPerIPU();

  Tensor Seqs = graph.addVariable(INT, {activeTiles, static_cast<size_t>(std::ceil(static_cast<double>(bufSize) / 4.0))}, "Seqs");

  // Metadata structure
  Tensor CompMeta = graph.addVariable(INT, {activeTiles, maxBatches * 4}, "CompMeta");

  auto [m, n] = similarityData.shape();

  poplar::Type sType;
  int workerMultiplier = 1;
  switch (vtype) {
    case VertexType::cpp: sType = INT; break;
    case VertexType::assembly: sType = FLOAT; break;
    case VertexType::multi: sType = INT; workerMultiplier = target.getNumWorkerContexts(); break;
    case VertexType::multiasm: sType = FLOAT; workerMultiplier = target.getNumWorkerContexts(); break;
    case VertexType::multistriped: sType = INT; break;  // all threads are going to work on the same problem
    case VertexType::multistripedasm: sType = FLOAT; break;  // all threads are going to work on the same problem
    case VertexType::stripedasm: sType = HALF; break;
  }

  TypeTraits traits = typeToTrait(sType);
  void* similarityBuffer;
  convertSimilarityMatrix(target, sType, similarityData, &similarityBuffer);
  Tensor similarity;
  if (vtype == VertexType::stripedasm) {
    similarity = graph.addConstant(sType, {m * n}, similarityBuffer, traits, false, "similarity");
  } else {
    similarity = graph.addConstant(sType, {m, n}, similarityBuffer, traits, false, "similarity");
  }
  free(similarityBuffer);

  Tensor Scores = graph.addVariable(INT, {activeTiles, maxBatches}, "Scores");
  Tensor ARanges = graph.addVariable(INT, {activeTiles, maxBatches}, "ARanges");
  Tensor BRanges = graph.addVariable(INT, {activeTiles, maxBatches}, "BRanges");

  graph.setTileMapping(similarity, 0);
  for (int i = 0; i < activeTiles; ++i) {
    int tileIndex = i % tileCount;
    graph.setTileMapping(Seqs[i], tileIndex);
    graph.setTileMapping(CompMeta[i], tileIndex);

    graph.setTileMapping(Scores[i], tileIndex);
    graph.setTileMapping(ARanges[i], tileIndex);
    graph.setTileMapping(BRanges[i], tileIndex);
  }

  OptionFlags streamOptions({/*{"bufferingDepth", "2"}, {"splitLimit", "0"}*/});
  // OptionFlags streamOptions({{"bufferingDepth", "2"}, {"splitLimit", "0"}});

  auto frontCs = graph.addComputeSet("SmithWaterman");
  for (int i = 0; i < activeTiles; ++i) {
    int tileIndex = i % tileCount;
    VertexRef vtx = graph.addVertex(frontCs, vertexTypeToString(vtype),
                                    {
                                        {"bufSize", bufSize},
                                        {"maxAB", maxAB},
                                        {"gapInit", gapInit},
                                        {"gapExt", gapExt},
                                        {"maxNPerTile", maxBatches},
                                        {"Seqs", Seqs[i]},
                                        {"Meta", CompMeta[i]},
                                        // {"Alen", Alens[i]},
                                        // {"Blen", Blens[i]},
                                        {"simMatrix", similarity},
                                        {"score", Scores[i]},
                                        {"ARange", ARanges[i]},
                                        {"BRange", BRanges[i]},
                                    });
    if (vtype == VertexType::stripedasm) {
      graph.connect(vtx["simWidth"], m);
      graph.setFieldSize(vtx["tS"], maxAB * workerMultiplier);
    }
    graph.setFieldSize(vtx["C"], maxAB * workerMultiplier);
    graph.setFieldSize(vtx["bG"], maxAB * workerMultiplier);
    if (vtype == VertexType::multistriped || vtype == VertexType::multistripedasm) {
      graph.setFieldSize(vtx["tS"], maxAB * target.getNumWorkerContexts());
      graph.setFieldSize(vtx["locks"], target.getNumWorkerContexts());
    }
    graph.setTileMapping(vtx, tileIndex);
    graph.setPerfEstimate(vtx, 1);
  }

  auto inputs_tensor = concat({Seqs.flatten(), CompMeta.flatten()});
  auto outputs_tensor = concat({Scores.flatten(), ARanges.flatten(), BRanges.flatten()});

  program::Sequence h2d_prog_concat;
  program::Sequence d2h_prog_concat;
  if (use_remote_buffer) {
    const char* units[] = {"B", "kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
    int b_size = inputs_tensor.numElements() * 4;
    PLOGI.printf("Use a RemoteBuffer of bytes %.2f", (double) b_size);
    int i = 0;
    for ( ;b_size > 1024; i++) {
        b_size /= 1024;
    }
   PLOGI.printf("Use a RemoteBuffer of banksize %.2f %s", (double) b_size, units[i]);
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
    auto host_stream_offset = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT, INT, 1);
    h2d_prog_concat.add(poplar::program::Copy(host_stream_offset, offset));
    h2d_prog_concat.add(program::PrintTensor(offset));
    auto remote_mem = graph.addRemoteBuffer(REMOTE_MEMORY, inputs_tensor.elementType(), inputs_tensor.numElements(), buf_rows, true, true);
    h2d_prog_concat.add(poplar::program::Copy(remote_mem, inputs_tensor, offset,"Copy Inputs from Remote->IPU"));
    // d2h_prog_concat.add(program::Execute(cs));

    auto device_stream_concat = graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL, INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements());
    d2h_prog_concat.add(poplar::program::Copy(outputs_tensor, device_stream_concat, "Copy Outputs from IPU->Host"));
  } else {
    auto host_stream_concat = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT, INT, inputs_tensor.numElements());
    auto device_stream_concat = graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL, INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements());
    h2d_prog_concat.add(poplar::program::Copy(host_stream_concat, inputs_tensor));
    d2h_prog_concat.add(poplar::program::Copy(outputs_tensor, device_stream_concat));
  }


  auto print_tensors_prog = program::Sequence({
    // program::PrintTensor("Alens", Alens),
    // program::PrintTensor("Blens", Blens),
    program::PrintTensor("CompMeta", CompMeta),
    // program::PrintTensor("Scores", Scores),
  });
#ifdef IPUMA_DEBUG
 program::Sequence main_prog;
  main_prog.add(program::Execute(frontCs));
  addCycleCount(graph, main_prog, CYCLE_COUNT_INNER);
#else
  auto main_prog = program::Execute(frontCs);
#endif
  prog.add(h2d_prog_concat);
  prog.add(main_prog);
  prog.add(d2h_prog_concat);

#ifdef IPUMA_DEBUG
  addCycleCount(graph, prog, CYCLE_COUNT_OUTER);
#endif
  return {prog, d2h_prog_concat, print_tensors_prog};
}

SWAlgorithm::SWAlgorithm(ipu::SWConfig config, IPUAlgoConfig algoconfig, int thread_id, bool useRemoteBuffer, size_t bufSize) : IPUAlgorithm(config), algoconfig(algoconfig), thread_id(thread_id), use_remote_buffer(useRemoteBuffer) {
  const auto totalComparisonsCount = algoconfig.getTotalNumberOfComparisons();
  buf_size = bufSize;
  buf_cap = buf_size;
  slots.resize(buf_size);
  std::fill(slots.begin(), slots.end(), false);

  scores.resize(totalComparisonsCount);
  a_range_result.resize(totalComparisonsCount);
  b_range_result.resize(totalComparisonsCount);

  Graph graph = createGraph();

  auto similarityMatrix = swatlib::selectMatrix(config.similarity, config.matchValue, config.mismatchValue, config.ambiguityValue);
  std::vector<program::Program> programs =
      buildGraph(graph, algoconfig.vtype, algoconfig.tilesUsed, algoconfig.maxAB, algoconfig.bufsize, algoconfig.maxBatches,
                 similarityMatrix, config.gapInit, config.gapExtend, use_remote_buffer, buf_size);

  createEngine(graph, programs);
}

SWAlgorithm::SWAlgorithm(ipu::SWConfig config, IPUAlgoConfig algoconfig) : SWAlgorithm::SWAlgorithm(config, algoconfig, 0) {}

BlockAlignmentResults SWAlgorithm::get_result() { return {scores, a_range_result, b_range_result}; }

void SWAlgorithm::checkSequenceSizes(const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B) {
  if (A.size() != B.size()) {
    PLOGW << "Mismatched size A " << A.size() << " != B " << B.size();
    exit(1);
  }
  if (A.size() > algoconfig.maxBatches * algoconfig.tilesUsed) {
    PLOGW << "A has more elements than the maxBatchsize";
    PLOGW << "A.size() = " << A.size();
    PLOGW << "max comparisons = " << algoconfig.tilesUsed << " * " << algoconfig.maxBatches;
    exit(1);
  }

  for (const auto& a : A) {
    if (a.size() > algoconfig.maxAB) {
      PLOGW << "Sequence size in a " << a.size() << " > " << algoconfig.maxAB;
      exit(1);
    }
  }
  for (const auto& b : B) {
    if (b.size() > algoconfig.maxAB) {
      PLOGW << "Sequence size in a " << b.size() << " > " << algoconfig.maxAB;
      exit(1);
    }
  }
}

std::vector<std::tuple<int, int>> SWAlgorithm::fillBuckets(const IPUAlgoConfig& algoconfig,const std::vector<std::string>& A, const std::vector<std::string>& B, int& err) {
  std::vector<std::tuple<int, int>> bucket_pairs;
  switch (algoconfig.fillAlgo) {
  case partition::Algorithm::fillFirst:
    err = partition::fillFirst(bucket_pairs, A, B, algoconfig.tilesUsed, algoconfig.bufsize, algoconfig.maxBatches);
    break;
  case partition::Algorithm::roundRobin:
    err = partition::roundRobin(bucket_pairs, A, B, algoconfig.tilesUsed, algoconfig.bufsize, algoconfig.maxBatches);
    break;
  case partition::Algorithm::greedy:
    err = partition::greedy(bucket_pairs, A, B, algoconfig.tilesUsed, algoconfig.bufsize, algoconfig.maxBatches);
    break;
  }
  return bucket_pairs;
}

size_t SWAlgorithm::getSeqsOffset(const IPUAlgoConfig& config) {
  return 0;
}

size_t SWAlgorithm::getMetaOffset(const IPUAlgoConfig& config) {
  return config.getSequenceBufferSize32b();
}

void SWAlgorithm::refetch() {
  engine->run(1);
}

void SWAlgorithm::upload(int32_t* inputs_begin, int32_t* inputs_end, slotToken slot) {
    PLOGD.printf("Slot is %d", slot);

    swatlib::TickTock rbt;
    rbt.tick();
    engine->copyToRemoteBuffer(inputs_begin, REMOTE_MEMORY, slot);
    release_slot(slot);
    rbt.tock();
    auto transferTime = rbt.duration<std::chrono::milliseconds>();
    size_t totalTransferSize = algoconfig.getInputBufferSize32b() * 4;
    auto transferBandwidth = (double) totalTransferSize / ((double) transferTime / 1000.0) / 1e6;
    PLOGI.printf("Transfer rate to remote buffer %.3f mb/s. %ld bytes in %.2fms", transferBandwidth, totalTransferSize, transferTime);
}

void SWAlgorithm::prepared_remote_compare(int32_t* inputs_begin, int32_t* inputs_end, int32_t* results_begin, int32_t* results_end, slotToken slot_token) {
  // We have to reconnect the streams to new memory locations as the destination will be in a shared memroy region.
  // engine->connectStream(HOST_STREAM_CONCAT, inputs_begin);
  engine->connectStream(STREAM_CONCAT_ALL, results_begin);

  std::vector<int> vv{slot_token};
  if (use_remote_buffer) {
    engine->connectStream(HOST_STREAM_CONCAT, vv.data());
  } else {
    engine->connectStream(HOST_STREAM_CONCAT, inputs_begin);
  }
  swatlib::TickTock t;
  t.tick();
  engine->run(0, "Execute Main");
  t.tock();
  release_slot(slot_token);
  PLOGD << "Total engine run time (in s): "  << static_cast<double>(t.duration<std::chrono::milliseconds>()) / 1000.0;

#ifdef IPUMA_DEBUG
  auto cyclesOuter = getTotalCycles(*engine, CYCLE_COUNT_OUTER);
  auto cyclesInner = getTotalCycles(*engine, CYCLE_COUNT_INNER);
  auto timeOuter = static_cast<double>(cyclesOuter) / getTarget().getTileClockFrequency();
  auto timeInner = static_cast<double>(cyclesInner) / getTarget().getTileClockFrequency();
  PLOGD << "Poplar cycle count: " << cyclesInner << "/" << cyclesOuter << " computed time (in s): " << timeInner << "/" << timeOuter;

  int32_t* meta_input = inputs_begin + getMetaOffset(algoconfig);

  // GCUPS computation
  uint64_t cellCount = 0;
  uint64_t dataCount = 0;
  for (size_t i = 0; i < algoconfig.getTotalNumberOfComparisons(); i++) {
    auto a_len = meta_input[4*i];
    auto b_len = meta_input[4*i+2];
    cellCount += a_len * b_len;
    dataCount += a_len + b_len;
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
  PLOGD << "Transfer time: " << transferTime << "s transfer ratio: " << transferInfoRatio << "% estimated bandwidth: " << transferBandwidth << "mb/s, per vertex: " << transferBandwidth / algoconfig.tilesUsed << "mb/s";
#endif
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
  prepare_remote(config, algoconfig, A, B, &*inputs.begin() + 2, &*inputs.end() - 2, mapping.data());
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

  // std::vector<int32_t> unord_scores(scores), unord_mismatches(mismatches), unord_a_range(a_range_result), unord_b_range(b_range_result);
  results[0] = 0xDEADBEEF;
  results[1] = 0xDEADBEEF;
  results[results_size + (1) + 1] = 0xFEEBDAED;
  results[results_size + (1) + 2] = 0xFEEBDAED;
  // prepared_remote_compare(a.data(), a_len.data(), b.data(), b_len.data(), unord_scores.data(), unord_mismatches.data(), unord_a_range.data(),
  //                         unord_b_range.data());
  slotToken slot = 0;
  if (use_remote_buffer) {
    slot = upload(&*inputs.begin() + 2, &*inputs.end() - 2);
  }
  tCompare.tick();
  prepared_remote_compare(&*inputs.begin() + 2, &*inputs.end() - 2, &*results.begin() + 2, &*results.end() - 2, slot);
  tCompare.tock();

  double prepare_time = static_cast<double>(tPrepare.duration<std::chrono::milliseconds>()) / 1e3;
  double compare_time = static_cast<double>(tCompare.duration<std::chrono::milliseconds>()) / 1e3;
  double total_time = prepare_time + compare_time;
  double prepare_perc = prepare_time / total_time * 100;
  PLOGD << "Total time: " << total_time << " prepare(" << prepare_time << ") compare(" << compare_time << ") preprocessing is " << prepare_perc << "% of total";

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

  // PLOGD << "Unordered results: " << swatlib::printVector(std::vector<int32_t>(results.begin() + 2, results.begin() + scores.size() + 2));

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
        PLOGW.printf("ERROR Thread %d received wrong data FIRST, try again data=%d, map_translate=%d\n", thread_id, scores[i], mapping[i] + 2);
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
    if ((double)sc/A.size() < 0.5) {
      PLOGW << "ERROR Too many scores are 0, retry number " << (nthTry - 1);
      refetch();
      goto retry;
    }
  }
}

std::string SWAlgorithm::printTensors() {
  std::stringstream graphstream;
  engine->setPrintTensorStream(graphstream);
  engine->run(2);
  return graphstream.str();
}

struct BucketFillInfo {
  size_t bucketContent;
  size_t bucketCmps;
  std::vector<int8_t> seqs;

  void resizeSeqs(int numSeqs) {
    seqs.resize(numSeqs, 0);
  }
};

std::vector<BucketMapping> SWAlgorithm::fillMNBuckets(const IPUAlgoConfig& algoconfig, const std::vector<std::string>& Seqs, const std::vector<int>& comparisons) {
  if (comparisons.size() & 0b1) {
    throw std::logic_error("Length of comparisons is not a multiple of 2, got " + std::to_string(comparisons.size()) + " instead");
  }
  int numComparisons = comparisons.size() >> 1;

  std::vector<BucketMapping> buckets(algoconfig.tilesUsed);
  std::vector<BucketFillInfo> bucketsInfo(algoconfig.tilesUsed);
  for (auto& b : bucketsInfo) {
    b.resizeSeqs(Seqs.size());
  }

  int bucketIndex = 0;
  bool incrementBucket = false;
  for (int n = 0; n < numComparisons; ++n) {
    int ai = comparisons[2*n];
    int bi = comparisons[2*n+1];

    int eff_asize, eff_bsize;
    bool mapped = false;
    for (int i = 0; i < buckets.size(); ++i) {
      int wrappedIndex = (bucketIndex + i) % buckets.size();
      auto& bInfo = bucketsInfo[wrappedIndex];
      auto& bBucket = buckets[wrappedIndex];

      int nextCmps, nextContentSize;
      nextCmps = bInfo.bucketCmps + 1;
      if (nextCmps > algoconfig.maxBatches) {
        continue;
      }
      // calculate needed space for sequences
      eff_asize = Seqs[ai].size();
      if (bInfo.seqs[ai] == 1) eff_asize = 0;
      eff_bsize = Seqs[bi].size();
      if (bInfo.seqs[bi] == 1) eff_bsize = 0;
      
      nextContentSize = bInfo.bucketContent + eff_asize + eff_bsize;
      if (nextContentSize > algoconfig.bufsize) {
        continue;
      }

      bInfo.bucketContent = nextContentSize;
      bInfo.bucketCmps = nextCmps;
      bInfo.seqs[ai] = 1;
      bInfo.seqs[bi] = 1;
      bBucket.comparisons.push_back({n, ai, bi});
      mapped = true;
      break;
    }
    if (!mapped) {
      throw std::runtime_error("Could not map comparison into a bucket");
    }
    if (eff_asize == 0 || eff_bsize == 0) {
      incrementBucket = false;
    }
    if (incrementBucket) {
      bucketIndex++;
      incrementBucket = false;
    } else {
      incrementBucket = true;
    }
  }
  return buckets;
}


void SWAlgorithm::transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end) {
  int maxComparisons = scores_end - scores_begin;
  size_t a_range_offset = maxComparisons;
  size_t b_range_offset = maxComparisons * 2;
  auto* results_score = results_begin;
  auto* results_arange = results_begin + a_range_offset;
  auto* results_brange = results_begin + b_range_offset;

  int numComparisons = mapping_end - mapping_begin;
  for (int i = 0; i < numComparisons; ++i) {
    auto mapped_i = mapping_begin[i];
    scores_begin[i] = results_score[mapped_i];
    arange_begin[i] = results_arange[mapped_i];
    brange_begin[i] = results_brange[mapped_i];
  }
}

void SWAlgorithm::compare_mn_local(const std::vector<std::string>& Seqs, const std::vector<int>& comparisons, bool errcheck) {
  size_t inputs_size = algoconfig.getInputBufferSize32b();
  std::vector<int32_t> inputs(inputs_size);
  size_t results_size = scores.size() + a_range_result.size() + b_range_result.size();
  std::vector<int32_t> results(results_size);

  auto numComparisons = comparisons.size() >> 1;
  auto cmpMapping = fillMNBuckets(algoconfig, Seqs, comparisons);
  auto mapping = fill_input_buffer(config, algoconfig, Seqs, cmpMapping, numComparisons, &*inputs.begin(), &*inputs.end());
  // std::cout << "Mapping: " << swatlib::printVector(mapping) << "\n";

  int slot = 0;
  if (use_remote_buffer) {
    slot = upload(&*inputs.begin(), &*inputs.end());
  }
  prepared_remote_compare(&*inputs.begin(), &*inputs.end(), &*results.begin(), &*results.end(), slot);

  transferResults(&*results.begin(), &*results.end(), &*mapping.begin(), &*mapping.end(), &*scores.begin(), &*scores.end(), &*a_range_result.begin(), &*a_range_result.end(), &*b_range_result.begin(), &*b_range_result.end());
}

struct BucketSequenceInfo {
  int bucketIndex;
  int bucketOffset;
};

BucketSequenceInfo& fillIfNotExists(std::unordered_map<int, BucketSequenceInfo>& bM, int sIndex, const std::string& data, int& bucketSeqN, int& bucketSeqOffset, int8_t* bSeq, swatlib::Encoding& encoder) {
  auto bit = bM.find(sIndex);
  // insert ai into bucket
  if (bit == bM.end()) {
    bM[sIndex] = {.bucketIndex = bucketSeqN, .bucketOffset = bucketSeqOffset};
    for (int k = 0; k < data.size(); ++k) {
      bSeq[bucketSeqOffset + k] = encoder.encode(data[k]);
    }
    bucketSeqN++;
    bucketSeqOffset += data.size();
  }
  return bM[sIndex];
}

std::vector<int> SWAlgorithm::fill_input_buffer(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& Seqs, const std::vector<BucketMapping>& comparisonMapping, int numComparisons, int32_t* inputs_begin, int32_t* inputs_end) {
  size_t input_elems = inputs_end - inputs_begin;
  memset(inputs_begin, 0, input_elems * sizeof(int32_t));

  int8_t* seqsInput = (int8_t*)inputs_begin + getSeqsOffset(algoconfig);
  int32_t* metaInput = inputs_begin + getMetaOffset(algoconfig);

  auto encoder = swatlib::getEncoder(swconfig.datatype);

  std::vector<int> mapping(numComparisons);

  for (int i = 0; i < comparisonMapping.size(); ++i) {
    auto& [comparisons] = comparisonMapping[i];
    std::unordered_map<int, BucketSequenceInfo> bucketSequences;
    int8_t* bSeq = seqsInput + i * algoconfig.getBufsize32b() * 4;
    int32_t* bMeta = metaInput + i * algoconfig.maxBatches * 4;

    int bucketSeqN = 0;
    int bucketSeqOffset = 0;

    for (int n = 0; n < comparisons.size(); ++n) {
      const auto [ci, ai, bi] = comparisons[n];
      auto& aInfo = fillIfNotExists(bucketSequences, ai, Seqs[ai], bucketSeqN, bucketSeqOffset, bSeq, encoder);
      auto& bInfo = fillIfNotExists(bucketSequences, bi, Seqs[bi], bucketSeqN, bucketSeqOffset, bSeq, encoder);
      bMeta[n*4 ] = Seqs[ai].size();
      bMeta[n*4+1] = aInfo.bucketOffset;
      bMeta[n*4+2] = Seqs[bi].size();
      bMeta[n*4+3] = bInfo.bucketOffset;

      if (aInfo.bucketOffset + Seqs[ai].size() > algoconfig.bufsize) {
        throw runtime_error("Offset larger than buffer size.");
      }
      if (bInfo.bucketOffset + Seqs[bi].size() > algoconfig.bufsize) {
        throw runtime_error("Offset larger than buffer size.");
      }

      // map results
      mapping[ci] = i * algoconfig.maxBatches + n;
    }
  }
  return mapping;
}

void SWAlgorithm::prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int32_t* inputs_begin, int32_t* inputs_end, int* seqMapping) {
  swatlib::TickTock preprocessTimer;
  std::vector<swatlib::TickTock> stageTimers(5);
  preprocessTimer.tick();
  stageTimers[0].tick();
  checkSequenceSizes(algoconfig, A, B);
  stageTimers[0].tock();

  stageTimers[1].tick();
  auto encoder = swatlib::getEncoder(swconfig.datatype);
  stageTimers[1].tock();

  stageTimers[2].tick();
  size_t input_elems = inputs_end - inputs_begin;
  memset(inputs_begin, 0, input_elems * sizeof(int32_t));
  stageTimers[2].tock();

  // #ifdef IPUMA_DEBUG
  // for (int32_t* it = inputs_begin; it != inputs_end; ++it) {
  //   if (*it != 0) {
  //     PLOGW << "Results are not zero";
  //   }
  // }
  // #endif

  stageTimers[3].tick();
  int errval = 0;
  auto mapping = fillBuckets(algoconfig, A, B, errval);

  if (errval) {
    PLOGW << "Bucket filling failed.";
    exit(1);
  }
  stageTimers[3].tock();

  stageTimers[4].tick();
  std::vector<std::tuple<int, int>> buckets(algoconfig.tilesUsed, {0, 0});

  const size_t seqs_offset = getSeqsOffset(algoconfig);
  const size_t meta_offset = getMetaOffset(algoconfig);

  int8_t* seqs = (int8_t*)inputs_begin + seqs_offset;
  int32_t* meta = inputs_begin + meta_offset;

  for (const auto [bucket, i] : mapping) {
    auto& [bN, bO] = buckets[bucket];
    const auto& ai = A[i];
    const auto& bi = B[i];
    const auto aSize = ai.size();
    const auto bSize = bi.size();

    const size_t offsetBuffer = bucket * algoconfig.getBufsize32b() * 4;
    const size_t offsetMeta = bucket * algoconfig.maxBatches * 4;
    auto* bucket_meta = meta + offsetMeta;
    auto* bucket_seq = seqs + offsetBuffer;

    if (bN*4+3 >= algoconfig.maxBatches * 4) {
      std::runtime_error("Bucket meta out of bounds.");
    }
    bucket_meta[bN*4] = aSize;
    bucket_meta[bN*4+1] = bO;
    if (bO + aSize > algoconfig.bufsize) {
      throw std::runtime_error("Sequence A too large for buffer.");
    }
    for (int j = 0; j < aSize; ++j) {
      bucket_seq[bO + j] = encoder.encode(ai[j]);
    }
    bO += aSize;

    bucket_meta[bN*4+2] = bSize;
    bucket_meta[bN*4+3] = bO;
    if (bO + bSize > algoconfig.bufsize) {
      throw std::runtime_error("Sequence B too large for buffer.");
    }
    for (int j = 0; j < bSize; ++j) {
      bucket_seq[bO + j] = encoder.encode(bi[j]);
    }
    bO += bSize;
    seqMapping[i] = bucket * algoconfig.maxBatches + bN;

    bN++;
  }
  stageTimers[4].tock();

  preprocessTimer.tock();
  PLOGD << "Total preprocessing time (in s): "  << static_cast<double>(preprocessTimer.duration<std::chrono::milliseconds>()) / 1000.0;

  PLOGD << "Stage timers:";
  for (int i = 0; i < stageTimers.size(); ++i)  {
    PLOGD << i << ": " << stageTimers[i].duration<std::chrono::milliseconds>();
  }

#ifdef IPUMA_DEBUG
  int emptyBuckets = 0;
  std::vector<int> bucketCmps;
  std::map<int, int> occurence;
  for (auto [n, bO] : buckets) {
    if (n == 0) emptyBuckets++;
    occurence[n]++;
    bucketCmps.push_back(n);
  }
  std::stringstream ss;
  ss << "Map[";
  for (auto [k ,v] : occurence) {
    ss << k << ": " << v << ",";
  }
  ss << "]";
  // SLOG(swatlib::printVector(bucketCmps), "\n");
  PLOGD << "Total number of buckets: " << buckets.size() << " empty buckets: " << emptyBuckets;
  PLOGD << "Bucket size occurence: " << ss.str();
#endif
  // SLOG("Inner comparison time: ", preprocessTimer.get_elapsed(), " engine run: ", engineTimer.get_elapsed(), "\n");
}
}  // namespace batchaffine
}  // namespace ipu