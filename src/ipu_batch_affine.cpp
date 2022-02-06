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

static const std::string CYCLE_COUNT_OUTER = "cycle-count-outer";
static const std::string CYCLE_COUNT_INNER = "cycle-count-inner";

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

int IPUAlgoConfig::getTotalNumberOfComparisons() const { return tilesUsed * maxBatches; }

int IPUAlgoConfig::getMetaBufferSize32b() const {return getTotalNumberOfComparisons() * 4;};
int IPUAlgoConfig::getLenBufferSize32b() const {return getTotalNumberOfComparisons() * 2;};

int IPUAlgoConfig::getSequenceBufferSize8b() const { return tilesUsed * bufsize; }
int IPUAlgoConfig::getSequenceBufferSize32b() const { return std::ceil(static_cast<double>(getSequenceBufferSize8b()) / 4.0); }

int IPUAlgoConfig::getInputBufferSize32b() const { return getSequenceBufferSize32b() * 2 + getMetaBufferSize32b(); }

std::string vertexTypeToString(VertexType v) { return typeString[static_cast<int>(v)]; }

/**
 * Streamable IPU graph for SW
 */
std::vector<program::Program> buildGraph(Graph& graph, VertexType vtype, unsigned long activeTiles, unsigned long maxAB,
                                         unsigned long bufSize, unsigned long maxBatches,
                                         const swatlib::Matrix<int8_t> similarityData, int gapInit, int gapExt) {
  program::Sequence prog;

  auto target = graph.getTarget();
  int tileCount = target.getTilesPerIPU();

  Tensor As = graph.addVariable(INT, {activeTiles, static_cast<size_t>(std::ceil(static_cast<double>(bufSize) / 4.0))}, "A");
  Tensor Bs = graph.addVariable(INT, {activeTiles, static_cast<size_t>(std::ceil(static_cast<double>(bufSize) / 4.0))}, "B");

  // Metadata structure
  Tensor CompMeta = graph.addVariable(INT, {activeTiles, maxBatches * 4}, "CompMeta");

  // Tensor Alens = graph.addVariable(INT, {activeTiles, maxBatches * 2}, "Alen");
  // Tensor Blens = graph.addVariable(INT, {activeTiles, maxBatches * 2}, "Blen");

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

  // Tensor similarity = graph.addConstant(SIGNED_CHAR, {m, n}, similarityData.data(), "similarity");
  Tensor Scores = graph.addVariable(INT, {activeTiles, maxBatches}, "Scores");
  Tensor ARanges = graph.addVariable(INT, {activeTiles, maxBatches}, "ARanges");
  Tensor BRanges = graph.addVariable(INT, {activeTiles, maxBatches}, "BRanges");

  graph.setTileMapping(similarity, 0);
  for (int i = 0; i < activeTiles; ++i) {
    int tileIndex = i % tileCount;
    graph.setTileMapping(As[i], tileIndex);
    graph.setTileMapping(Bs[i], tileIndex);
    graph.setTileMapping(CompMeta[i], tileIndex);
    // graph.setTileMapping(Alens[i], tileIndex);
    // graph.setTileMapping(Blens[i], tileIndex);

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
                                        {"A", As[i]},
                                        {"B", Bs[i]},
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

  auto inputs_tensor = concat({As.flatten(), Bs.flatten(), CompMeta.flatten()});
  auto outputs_tensor = concat({Scores.flatten(), ARanges.flatten(), BRanges.flatten()});

  auto host_stream_concat = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT, INT, inputs_tensor.numElements());
  auto device_stream_concat = graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL, INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements());
  auto h2d_prog_concat = program::Sequence({poplar::program::Copy(host_stream_concat, inputs_tensor)});
  auto d2h_prog_concat = program::Sequence({poplar::program::Copy(outputs_tensor, device_stream_concat)});

  auto print_tensors_prog = program::Sequence({
    // program::PrintTensor("Alens", Alens),
    // program::PrintTensor("Blens", Blens),
    program::PrintTensor("CompMeta", CompMeta),
    // program::PrintTensor("Scores", Scores),
  });
#ifdef IPUMA_DEBUG
 program::Sequence  main_prog;
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

SWAlgorithm::SWAlgorithm(ipu::SWConfig config, IPUAlgoConfig algoconfig, int thread_id) : IPUAlgorithm(config), algoconfig(algoconfig), thread_id(thread_id) {
  const auto totalComparisonsCount = algoconfig.getTotalNumberOfComparisons();

  scores.resize(totalComparisonsCount);
  a_range_result.resize(totalComparisonsCount);
  b_range_result.resize(totalComparisonsCount);

  Graph graph = createGraph();

  auto similarityMatrix = swatlib::selectMatrix(config.similarity, config.matchValue, config.mismatchValue, config.ambiguityValue);
  std::vector<program::Program> programs =
      buildGraph(graph, algoconfig.vtype, algoconfig.tilesUsed, algoconfig.maxAB, algoconfig.bufsize, algoconfig.maxBatches,
                 similarityMatrix, config.gapInit, config.gapExtend);

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

std::vector<std::tuple<int, int>> SWAlgorithm::fillBuckets(const std::vector<std::string>& A, const std::vector<std::string>& B, int& err) {
  return fillBuckets(algoconfig, A, B, err);
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

size_t SWAlgorithm::getAOffset(const IPUAlgoConfig& config) {
  return 0;
}

size_t SWAlgorithm::getBOffset(const IPUAlgoConfig& config) {
  return config.getSequenceBufferSize32b();
}

size_t SWAlgorithm::getMetaOffset(const IPUAlgoConfig& config) {
  return config.getSequenceBufferSize32b() * 2;
}

void SWAlgorithm::refetch() {
  engine->run(1);
}

void SWAlgorithm::prepared_remote_compare(int32_t* inputs_begin, int32_t* inputs_end, int32_t* results_begin, int32_t* results_end) {
  // We have to reconnect the streams to new memory locations as the destination will be in a shared memroy region.
  engine->connectStream(HOST_STREAM_CONCAT, inputs_begin);
  engine->connectStream(STREAM_CONCAT_ALL, results_begin);

  swatlib::TickTock t;
  t.tick();
  engine->run(0);
  t.tock();

#ifdef IPUMA_DEBUG
  auto cyclesOuter = getTotalCycles(*engine, CYCLE_COUNT_OUTER);
  auto cyclesInner = getTotalCycles(*engine, CYCLE_COUNT_INNER);
  auto timeOuter = static_cast<double>(cyclesOuter) / getTarget().getTileClockFrequency();
  auto timeInner = static_cast<double>(cyclesInner) / getTarget().getTileClockFrequency();
  PLOGD << "Poplar cycle count: " << cyclesInner << "/" << cyclesOuter << " computed time (in s): " << timeInner << "/" << timeOuter;

  // size_t alen_offset = getAlenOffset(algoconfig);
  // size_t blen_offset = getBlenOffset(algoconfig);

  int32_t* meta_input = inputs_begin + getMetaOffset(algoconfig);

  // int32_t* a_len = inputs_begin + alen_offset;
  // int32_t* b_len = inputs_begin + blen_offset;

  // GCUPS computation
  // auto cellCount = getCellCount(A, B);
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
  std::vector<int> mapping;
  size_t inputs_size = algoconfig.getInputBufferSize32b();
  std::vector<int32_t> inputs(inputs_size + 4);

  size_t results_size = scores.size() + a_range_result.size() + b_range_result.size();
  std::vector<int32_t> results(results_size + 4);

  inputs[0] = 0xDEADBEEF;
  inputs[1] = 0xDEADBEEF;
  inputs[inputs_size + (1) + 1] = 0xFEEBDAED;
  inputs[inputs_size + (1) + 2] = 0xFEEBDAED;
  prepare_remote(config, algoconfig, A, B, &*inputs.begin() + 2, &*inputs.end() - 2, mapping);

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
  prepared_remote_compare(&*inputs.begin() + 2, &*inputs.end() - 2, &*results.begin() + 2, &*results.end() - 2);

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
  // return;
// retry:
//   prepared_remote_compare(a.data(), a_len.data(), b.data(), b_len.data(), unord_scores.data(), unord_mismatches.data(), unord_a_range.data(),
//                           unord_b_range.data());
//   // reorder results based on mapping
//   for (int i = 0; i < mapping.size(); ++i) {
//     scores[i] = unord_scores[mapping[i]];
//     if (scores[i] >= KLIGN_IPU_MAXAB_SIZE) {
//       printf("Thread %d received wrong data AGAIN, FATAL data=%d, map_translate=%d\n", -1, scores[i], mapping[i]);
//       exit(1);
//     }
//     mismatches[i] = unord_mismatches[mapping[i]];
//     a_range_result[i] = unord_a_range[mapping[i]];
//     b_range_result[i] = unord_b_range[mapping[i]];
//   }
}

std::string SWAlgorithm::printTensors() {
  std::stringstream graphstream;
  engine->setPrintTensorStream(graphstream);
  engine->run(2);
  return graphstream.str();
}

void SWAlgorithm::prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int32_t* inputs_begin, int32_t* inputs_end, std::vector<int>& seqMapping) {
  swatlib::TickTock preprocessTimer;
  preprocessTimer.tick();
  checkSequenceSizes(algoconfig, A, B);

  auto encoder = swatlib::getEncoder(swconfig.datatype);
  auto vA = encoder.encode(A);
  auto vB = encoder.encode(B);

  size_t input_elems = inputs_end - inputs_begin;
  memset(inputs_begin, 0, input_elems * sizeof(int32_t));

  #ifdef IPUMA_DEBUG
  for (int32_t* it = inputs_begin; it != inputs_end; ++it) {
    if (*it != 0) {
      PLOGW << "Results are not zero";
    }
  }
  #endif

  int errval = 0;
  auto mapping = fillBuckets(algoconfig, A, B, errval);

  if (errval) {
    PLOGW << "Bucket filling failed.";
    exit(1);
  }

  std::vector<std::tuple<int, int, int>> buckets(algoconfig.tilesUsed, {0, 0, 0});

  seqMapping = std::vector<int>(A.size(), 0);
  size_t a_offset = getAOffset(algoconfig);
  size_t b_offset = getBOffset(algoconfig);
  size_t meta_offset = getMetaOffset(algoconfig);
  // size_t alen_offset = getAlenOffset(algoconfig);
  // size_t blen_offset = getBlenOffset(algoconfig);

  int8_t* a = (int8_t*)inputs_begin;
  int8_t* b = (int8_t*)(inputs_begin + b_offset);
  int32_t* meta = inputs_begin + meta_offset;
  // int32_t* a_len = inputs_begin + alen_offset;
  // int32_t* b_len = inputs_begin + blen_offset;

  for (const auto [bucket, i] : mapping) {
    auto& [bN, bA, bB] = buckets[bucket];
    auto aSize = vA[i].size();
    auto bSize = vB[i].size();

    size_t offsetBuffer = bucket * algoconfig.bufsize;
    size_t offsetMeta = bucket * algoconfig.maxBatches * 4;
    auto* bucket_meta = meta + offsetMeta;

    bucket_meta[bN*4  ] = aSize;
    bucket_meta[bN*4+1] = bA;
    bucket_meta[bN*4+2] = bSize;
    bucket_meta[bN*4+3] = bB;
    seqMapping[i] = bucket * algoconfig.maxBatches + bN;

    memcpy(a + offsetBuffer + bA, vA[i].data(), aSize);
    memcpy(b + offsetBuffer + bB, vB[i].data(), bSize);

    bN++;
    bA += aSize;
    bB += bSize;
  }

  preprocessTimer.tock();

#ifdef IPUMA_DEBUG
  int emptyBuckets = 0;
  std::vector<int> bucketCmps;
  std::map<int, int> occurence;
  for (auto [n, bA, bB] : buckets) {
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