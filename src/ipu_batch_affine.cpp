#include <string>
#include <cmath>
#include <iostream>

#include <plog/Log.h>

#include <poplar/Graph.hpp>
#include <poputil/TileMapping.hpp>
#include <popops/Zero.hpp>

#include "swatlib/swatlib.h"
#include "ipu_batch_affine.h"

namespace ipu {
namespace batchaffine {

static const std::string STREAM_CONCAT_ALL = "concat-read-all";
static const std::string HOST_STREAM_CONCAT = "host-stream-concat";
static const std::string REMOTE_MEMORY = "inputs-remotebuffer";

static const std::string CYCLE_COUNT_OUTER = "cycle-count-outer";
static const std::string CYCLE_COUNT_INNER = "cycle-count-inner";

/**
 * Streamable IPU graph for SW
 */
std::vector<program::Program> buildGraph(Graph& graph, VertexType vtype, unsigned long activeTiles, unsigned long maxAB,
                                         unsigned long bufSize, unsigned long maxBatches,
                                         const swatlib::Matrix<int8_t> similarityData, int gapInit, int gapExt,
                                         bool use_remote_buffer, int buf_rows, bool forward_only) {
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

  OptionFlags streamOptions({});

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
    graph.setFieldSize(vtx["C"], maxAB * workerMultiplier);
    graph.setFieldSize(vtx["bG"], maxAB * workerMultiplier);
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
    auto host_stream_offset = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT, INT, 1);
    h2d_prog_concat.add(poplar::program::Copy(host_stream_offset, offset));
    auto remote_mem =
        graph.addRemoteBuffer(REMOTE_MEMORY, inputs_tensor.elementType(), inputs_tensor.numElements(), buf_rows, true, true);
    h2d_prog_concat.add(poplar::program::Copy(remote_mem, inputs_tensor, offset, "Copy Inputs from Remote->IPU"));
    // d2h_prog_concat.add(program::Execute(cs));

    auto device_stream_concat =
        graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL, INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements());
    d2h_prog_concat.add(poplar::program::Copy(outputs_tensor, device_stream_concat, "Copy Outputs from IPU->Host"));
  } else {
    auto host_stream_concat = graph.addHostToDeviceFIFO(HOST_STREAM_CONCAT, INT, inputs_tensor.numElements());
    auto device_stream_concat =
        graph.addDeviceToHostFIFO(STREAM_CONCAT_ALL, INT, Scores.numElements() + ARanges.numElements() + BRanges.numElements());
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

SWAlgorithm::SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, bool useRemoteBuffer, size_t slotCap)
    : IPUAlgorithm(config)
    , algoconfig(algoconfig)
    , thread_id(thread_id)
    , use_remote_buffer(useRemoteBuffer) {
  const auto totalComparisonsCount = algoconfig.getTotalNumberOfComparisons();
  slot_avail = slotCap;
  slots.resize(slot_avail);
  std::fill(slots.begin(), slots.end(), false);

  scores.resize(totalComparisonsCount);
  a_range_result.resize(totalComparisonsCount);
  b_range_result.resize(totalComparisonsCount);

  Graph graph = createGraph();

  auto similarityMatrix = swatlib::selectMatrix(config.similarity, config.matchValue, config.mismatchValue, config.ambiguityValue);
  std::vector<program::Program> programs =
      buildGraph(graph, algoconfig.vtype, algoconfig.tilesUsed, algoconfig.maxAB, algoconfig.bufsize, algoconfig.maxBatches,
                 similarityMatrix, config.gapInit, config.gapExtend, use_remote_buffer, slot_avail, algoconfig.forwardOnly);

  createEngine(graph, programs);
}

SWAlgorithm::SWAlgorithm(ipu::SWConfig config, IPUAlgoConfig algoconfig)
    : SWAlgorithm::SWAlgorithm(config, algoconfig, 0) {}

BlockAlignmentResults SWAlgorithm::get_result() { return {scores, a_range_result, b_range_result}; }

void SWAlgorithm::checkSequenceSizes(const IPUAlgoConfig& algoconfig, const std::vector<int>& SeqSizes) {
  if (SeqSizes.size() > algoconfig.maxBatches * algoconfig.tilesUsed) {
    PLOGW << "Sequence has more elements than the maxBatchsize";
    PLOGW << "Seq.size() = " << SeqSizes.size();
    PLOGW << "max comparisons = " << algoconfig.tilesUsed << " * " << algoconfig.maxBatches;
    exit(1);
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
  std::transform(Seqs.begin(), Seqs.end(), std::back_inserter(seqSizes), [](const auto& s) {return s.size();});
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

size_t SWAlgorithm::getMetaOffset(const IPUAlgoConfig& config) { return config.getTotalBufsize32b(); }

void SWAlgorithm::refetch() { engine->run(1); }

void SWAlgorithm::upload(int32_t* inputs_begin, int32_t* inputs_end, slotToken slot) {
  PLOGD.printf("Slot is %d", slot);
  swatlib::TickTock rbt;
  rbt.tick();
  engine->copyToRemoteBuffer(inputs_begin, REMOTE_MEMORY, slot);
  rbt.tock();
  auto transferTime = rbt.duration<std::chrono::milliseconds>();
  size_t totalTransferSize = algoconfig.getInputBufferSize32b() * 4;
  auto transferBandwidth = (double)totalTransferSize / ((double)transferTime / 1000.0) / 1e6;
  PLOGI.printf("Transfer rate to remote buffer %.3f mb/s. %ld bytes in %.2fms", transferBandwidth, totalTransferSize, transferTime);
}

void SWAlgorithm::prepared_remote_compare(int32_t* inputs_begin, int32_t* inputs_end, int32_t* results_begin, int32_t* results_end,
                                          slotToken slot_token) {
  // We have to reconnect the streams to new memory locations as the destination will be in a shared memroy region.
  // engine->connectStream(HOST_STREAM_CONCAT, inputs_begin);
  engine->connectStream(STREAM_CONCAT_ALL, results_begin);

  std::vector<int> vv{slot_token};
  if (use_remote_buffer) {
    PLOGD.printf("Use remote buffer with slot %d", vv[0]);
    engine->connectStream(HOST_STREAM_CONCAT, vv.data());
  } else {
    engine->connectStream(HOST_STREAM_CONCAT, inputs_begin);
  }
  swatlib::TickTock t;
  t.tick();
  engine->run(0, "Execute Main");
  t.tock();
  release_slot(slot_token);
  PLOGD << "Total engine run time (in s): " << static_cast<double>(t.duration<std::chrono::milliseconds>()) / 1000.0;

#ifdef IPUMA_DEBUG
  auto cyclesOuter = getTotalCycles(*engine, CYCLE_COUNT_OUTER);
  auto cyclesInner = getTotalCycles(*engine, CYCLE_COUNT_INNER);
  auto timeOuter = static_cast<double>(cyclesOuter) / getTarget().getTileClockFrequency();
  auto timeInner = static_cast<double>(cyclesInner) / getTarget().getTileClockFrequency();
  PLOGD << "Poplar cycle count: " << cyclesInner << "/" << cyclesOuter << " computed time (in s): " << timeInner << "/"
        << timeOuter;

  int32_t* meta_input = inputs_begin + getMetaOffset(algoconfig);

  // GCUPS computation
  uint64_t cellCount = 0;
  uint64_t dataCount = 0;
  for (size_t i = 0; i < algoconfig.getTotalNumberOfComparisons(); i++) {
    auto a_len = meta_input[4 * i];
    auto b_len = meta_input[4 * i + 2];
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
  PLOGD << "Transfer time: " << transferTime << "s estimated bandwidth: " << transferBandwidth
        << "mb/s, per vertex: " << transferBandwidth / algoconfig.tilesUsed << "mb/s";
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
  std::stringstream graphstream;
  engine->setPrintTensorStream(graphstream);
  engine->run(2);
  return graphstream.str();
}

void SWAlgorithm::transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end) {
  int numComparisons = mapping_end - mapping_begin;
  transferResults(results_begin, results_end, mapping_begin, mapping_end, scores_begin, scores_end, arange_begin, arange_end, brange_begin, brange_end);
}
void SWAlgorithm::transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end, int numComparisons) {
  size_t results_size = results_end - results_begin;
  size_t result_part_size = results_size / 3;
  if (results_size % 3 != 0)  {
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

#ifdef IPUMA_DEBUG
  int emptyBuckets = 0;
  int maxBucket = 0;
  long dataCount = 0;
  std::vector<int> bucketCmps;
  std::map<int, int> occurence;
  for (const auto& bucket : map.buckets) {
    if (bucket.cmps.size() == 0) emptyBuckets++;
    occurence[bucket.cmps.size()]++;
    bucketCmps.push_back(bucket.cmps.size());
    maxBucket = std::max(maxBucket, bucket.seqSize);
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
  double bucketPerc = static_cast<double>(maxBucket) / static_cast<double>(algoconfig.bufsize) * 100.0;
  PLOGD << "Max bucket: " << maxBucket << "/" << algoconfig.bufsize << " (" << bucketPerc << "%)\n";
  double totalTransferSize = algoconfig.getInputBufferSize32b() * 4;
  auto transferInfoRatio = static_cast<double>(dataCount) / totalTransferSize * 100;
  PLOGD << "Transfer info/total: " << dataCount << "/" << totalTransferSize << " (" << transferInfoRatio << "%)\n";
#endif

  int slot = 0;
  if (use_remote_buffer) {
    slot = upload(&*inputs.begin(), &*inputs.end());
  }
  prepared_remote_compare(&*inputs.begin(), &*inputs.end(), &*results.begin(), &*results.end(), slot);

  transferResults(&*results.begin(), &*results.end(), &*mapping.begin(), &*mapping.end(), &*scores.begin(), &*scores.end(),
                  &*a_range_result.begin(), &*a_range_result.end(), &*b_range_result.begin(), &*b_range_result.end());
}

void SWAlgorithm::fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& Seqs, const Comparisons& Cmps, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping) {
  auto encodeTable = swatlib::getEncoder(dtype).getCodeTable();
  const size_t seqs_offset = getSeqsOffset(algoconfig);
  const size_t meta_offset = getMetaOffset(algoconfig);

  int8_t* seqs = (int8_t*)inputs_begin + seqs_offset;
  int32_t* meta = inputs_begin + meta_offset;

  for (const auto& bucketMapping : map.buckets) {
    const size_t offsetBuffer = bucketMapping.bucketIndex * algoconfig.getBufsize32b() * 4;
    const size_t offsetMeta = bucketMapping.bucketIndex * algoconfig.maxBatches * 4;
    auto* bucket_meta = meta + offsetMeta;
    auto* bucket_seq = seqs + offsetBuffer;

    for (const auto& sm : bucketMapping.seqs) {
      const char* seq;
      int seqSize;
      switch(sm.origin) {
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
      bucket_meta[i*4  ] = cmpMapping.sizeA;
      bucket_meta[i*4+1] = cmpMapping.offsetA;
      bucket_meta[i*4+2] = cmpMapping.sizeB;
      bucket_meta[i*4+3] = cmpMapping.offsetB;

      mapping[cmpMapping.comparisonIndex] = map.cmpCapacity * bucketMapping.bucketIndex + i;
    }
  }
}

void SWAlgorithm::fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& A, const RawSequences& B, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping) {
  size_t input_elems = inputs_end - inputs_begin;
  memset(inputs_begin, 0, input_elems * sizeof(int32_t));

  const auto encodeTable = swatlib::getEncoder(dtype).getCodeTable();
  const size_t seqs_offset = getSeqsOffset(algoconfig);
  const size_t meta_offset = getMetaOffset(algoconfig);

  int8_t* seqs = (int8_t*)inputs_begin + seqs_offset;
  int32_t* meta = inputs_begin + meta_offset;

  for (const auto& bucketMapping : map.buckets) {
    const size_t offsetBuffer = bucketMapping.bucketIndex * algoconfig.getBufsize32b() * 4;
    const size_t offsetMeta = bucketMapping.bucketIndex * algoconfig.maxBatches * 4;
    auto* bucket_meta = meta + offsetMeta;
    auto* bucket_seq = seqs + offsetBuffer;

    for (const auto& sm : bucketMapping.seqs) {
      const char* seq;
      int seqSize;
      switch(sm.origin) {
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
      for (int j = 0; j < seqSize; ++j) {
        bucket_seq[sm.offset + j] = encodeTable[seq[j]];
      }
    }
    for (int i = 0; i < bucketMapping.cmps.size(); ++i) {
      const auto& cmpMapping = bucketMapping.cmps[i];
      bucket_meta[i*4  ] = cmpMapping.sizeA;
      bucket_meta[i*4+1] = cmpMapping.offsetA;
      bucket_meta[i*4+2] = cmpMapping.sizeB;
      bucket_meta[i*4+3] = cmpMapping.offsetB;

      mapping[cmpMapping.comparisonIndex] = map.cmpCapacity * bucketMapping.bucketIndex + i;
    }
  }
}

void SWAlgorithm::prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int32_t* inputs_begin, int32_t* inputs_end, int* seqMapping) {
  swatlib::TickTock preprocessTimer;
  std::vector<swatlib::TickTock> stageTimers(3);
  preprocessTimer.tick();
  checkSequenceSizes(algoconfig, A);
  checkSequenceSizes(algoconfig, B);

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

#ifdef IPUMA_DEBUG
  int emptyBuckets = 0;
  int maxBucket = 0;
  long dataCount = 0;
  std::vector<int> bucketCmps;
  std::map<int, int> occurence;
  for (const auto& bucket : map.buckets) {
    if (bucket.cmps.size() == 0) emptyBuckets++;
    occurence[bucket.cmps.size()]++;
    bucketCmps.push_back(bucket.cmps.size());
    maxBucket = std::max(maxBucket, bucket.seqSize);
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
  double bucketPerc = static_cast<double>(maxBucket) / static_cast<double>(algoconfig.bufsize) * 100.0;
  PLOGD << "Max bucket: " << maxBucket << "/" << algoconfig.bufsize << " (" << bucketPerc << "%)\n";
  double totalTransferSize = algoconfig.getInputBufferSize32b() * 4;
  auto transferInfoRatio = static_cast<double>(dataCount) / totalTransferSize * 100;
  PLOGD << "Transfer info/total: " << dataCount << "/" << totalTransferSize << " (" << transferInfoRatio << "%)\n";
#endif
}
slotToken SWAlgorithm::queue_slot() {
  assert(buf_has_capacity());
  int s = -1;
  for (size_t i = 0; i < slots.size(); i++) {
    if (!slots[i]) {
      s = i;
      break;
    }
  }
  assert(s != -1);
  slots[s] = true;
  slot_avail--;
  last_slot = s;
  return s;
}
void SWAlgorithm::release_slot(slotToken i) {
  assert(slots[i] == true);
  slots[i] = false;
  slot_avail++;
}

bool SWAlgorithm::slot_available() { return slot_avail > 0; }

slotToken SWAlgorithm::upload(int32_t* inputs_begin, int32_t* inputs_end) {
  int slot = queue_slot();
  upload(inputs_begin, inputs_end, slot);
  return slot;
}

}  // namespace batchaffine
}  // namespace ipu