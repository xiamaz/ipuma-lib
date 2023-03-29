#include "batching.h"
#include <plog/Log.h>
#include <omp.h>

#include "partition.h"

namespace ipu {

int32_t getComparisonIndex(OriginIndex oi) {
  return oi / NSEEDS;
}

int32_t getSeedIndex(OriginIndex oi) {
  return oi % NSEEDS;
}

std::tuple<int32_t, int32_t> unpackOriginIndex(OriginIndex oi) {
  return {getComparisonIndex(oi), getSeedIndex(oi)};
}

OriginIndex packOriginIndex(int32_t comparisonIndex, int32_t seedIndex) {
  return comparisonIndex * NSEEDS + seedIndex;
}

std::string Batch::toString() const {
  std::stringstream ss;
  ss << "Batch[" << this->numComparisons << " size: " << this->dataCount << "]";
  return ss.str();
}

int8_t* Batch::getSequenceBuffer() {
  return (int8_t*) inputs.data();
}

int8_t* Batch::getMetaBuffer() {
  return (int8_t*) (inputs.data() + metaOffset);
}

Batch::Batch() {
}

Batch::Batch(IPUAlgoConfig config) {
  initialize(config);
}
void Batch::clear() {
  inputs.resize(0);
  results.resize(0);
  origin_comparison_index.resize(0);
  metaOffset = 0;
  maxComparisons = 0;
  numComparisons = 0;
  cellCount = 0;
  dataCount = 0;
}

void Batch::initialize(IPUAlgoConfig config) {
  inputs.resize(config.getInputBufferSize32b(), 0);
  results.resize(config.getTotalNumberOfComparisons() * 3, 0);
  origin_comparison_index.resize(config.getTotalNumberOfComparisons(), -1 * NSEEDS);
  metaOffset = config.getOffsetMetadata();
  maxComparisons = config.getTotalNumberOfComparisons();
  numComparisons = 0;
  cellCount = 0;
  dataCount = 0;
}

BlockAlignmentResults Batch::get_result() {
  std::vector<int32_t> scores(maxComparisons);
  std::vector<uint32_t> a_range_result(maxComparisons);
  std::vector<uint32_t> b_range_result(maxComparisons);

  // TODO this should not be hardcoded here but passed as a config to make layout changes easier.
  int aOffset = results.size() / 3;
  int bOffset = aOffset * 2;

  // mapping is now responsibility of the caller
  for (int i = 0; i < maxComparisons; ++i) {
    int aindex = i + aOffset;
    int bindex = i + bOffset;
    scores[i] = results[i];
    a_range_result[i] = results[aOffset + i];
    b_range_result[i] = results[bOffset + i];
  }

  return {scores, a_range_result, b_range_result};
}

Batch create_batch(const partition::BatchMapping& map, const RawSequences& seqs, const IPUAlgoConfig& algoconfig, const SWConfig& config) {
  Batch batch;
  batch.initialize(algoconfig);

  auto encodeTable = swatlib::getEncoder(config.datatype).getCodeTable();
  int8_t* seqInput = batch.getSequenceBuffer();
  int8_t* metaInput = batch.getMetaBuffer();

  auto& cellCount = batch.cellCount;
  auto& dataCount = batch.dataCount;
  auto& comparisonCount = batch.numComparisons;

  for (int bi = 0; bi < map.buckets.size(); bi++) {
    const auto& bucketMapping = map.buckets[bi];
    const size_t offsetSequence = bucketMapping.bucketIndex * algoconfig.getBufsize32b() * 4;

    auto* bucketSeq = seqInput + offsetSequence;
    for (int si = 0; si< bucketMapping.seqs.size(); si++) {
      const auto& sequenceMapping = bucketMapping.seqs[si];
      const char *seq = seqs[sequenceMapping.index].data();
      size_t seqSize = seqs[sequenceMapping.index].size();
      #pragma omp simd
      for (int j = 0; j < seqSize; ++j) {
        bucketSeq[sequenceMapping.offset + j] = encodeTable[seq[j]];
      }
      dataCount += seqSize;
    }

    for (int i = 0; i < bucketMapping.cmps.size(); ++i) {
      const auto& comparisonMapping = bucketMapping.cmps[i];
      const auto& cmp = comparisonMapping.comparison;
      if (isSeeded) {
        auto* bucketMeta = (XDropMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
        // PLOGE << "Num Comparisons " << bucketMapping.cmps.size() << " " << NSEEDS;
        for (int j = 0; j < NSEEDS; ++j) {
          bucketMeta[i * NSEEDS + j] = {
            comparisonMapping.createMeta(),
            .seedAStartPos = static_cast<int32_t>(cmp.seeds[j].seedAStartPos),
            .seedBStartPos = static_cast<int32_t>(cmp.seeds[j].seedBStartPos),
          };
          // PLOGE << "Seeds " << bucketMeta[i * NSEEDS + j].seedAStartPos << " " <<  bucketMeta[i * NSEEDS + j].seedBStartPos;
          // PLOGE << "Lengths " << bucketMeta[i * NSEEDS + j].sizeA << " " <<  bucketMeta[i * NSEEDS + j].sizeB;
          auto validCmp = bucketMeta[i * NSEEDS + j].seedAStartPos != -1;

          batch.origin_comparison_index[algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex + i * NSEEDS + j] = packOriginIndex(cmp.originalComparisonIndex, j);
          cellCount += cmp.sizeA * cmp.sizeB * validCmp;
        }
      } else {
        auto* bucketMeta = (SWMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
        bucketMeta[i] = comparisonMapping.createMeta();
        batch.origin_comparison_index[algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex + i] = packOriginIndex(cmp.originalComparisonIndex, 0);
        cellCount += cmp.sizeA * cmp.sizeB;
      }

      comparisonCount++;
    }
  }
  return std::move(batch);
}

template<typename C>
std::vector<Batch> create_batches(const RawSequences& seqs, std::vector<C>& cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config) {
  std::vector<swatlib::TickTock> stageTimers(3);
  stageTimers[0].tick();
  stageTimers[1].tick();
  auto mappings = partition::mapBatches(algoconfig, seqs, cmps);
  stageTimers[1].tock();

  std::vector<Batch> batches(mappings.size());
  const auto inputBufferSize = algoconfig.getInputBufferSize32b();
  auto encodeTable = swatlib::getEncoder(config.datatype).getCodeTable();
  const bool isSeeded = algoconfig.hasSeeds();

  stageTimers[2].tick();
  size_t vertexMetaSize = algoconfig.maxComparisonsPerVertex * sizeof(SWMeta);
  if (isSeeded) {
    vertexMetaSize = algoconfig.maxComparisonsPerVertex * sizeof(XDropMeta);
  }

  swatlib::TickTock seqT;
  swatlib::TickTock cmpT;
  auto cmps_total = 0;

  #pragma omp parallel for schedule(static) num_threads(48)
  for (int mi = 0; mi < mappings.size(); ++mi) {
    const auto& map = mappings[mi];
    Batch& batch = batches[mi];
    batch.initialize(algoconfig);

    int8_t* seqInput = batch.getSequenceBuffer();
    int8_t* metaInput = batch.getMetaBuffer();
    auto& cellCount = batch.cellCount;
    auto& dataCount = batch.dataCount;
    auto& comparisonCount = batch.numComparisons;

    for (int bi = 0; bi < map.buckets.size(); bi++) {
      const auto& bucketMapping = map.buckets[bi];
      const size_t offsetSequence = bucketMapping.bucketIndex * algoconfig.getBufsize32b() * 4;

      auto* bucketSeq = seqInput + offsetSequence;
      seqT.tick();
      for (int si = 0; si< bucketMapping.seqs.size(); si++) {
        const auto& sequenceMapping = bucketMapping.seqs[si];
        const char *seq = seqs[sequenceMapping.index].data();
        size_t seqSize = seqs[sequenceMapping.index].size();
        #pragma omp simd
        for (int j = 0; j < seqSize; ++j) {
          bucketSeq[sequenceMapping.offset + j] = encodeTable[seq[j]];
        }
        dataCount += seqSize;
      }
      seqT.tock();

      cmpT.tick();
      for (int i = 0; i < bucketMapping.cmps.size(); ++i) {
        const auto& comparisonMapping = bucketMapping.cmps[i];
        const auto& cmp = comparisonMapping.comparison;
        if (isSeeded) {
          auto* bucketMeta = (XDropMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
          // PLOGE << "Num Comparisons " << bucketMapping.cmps.size() << " " << NSEEDS;
          for (int j = 0; j < NSEEDS; ++j) {
            bucketMeta[i * NSEEDS + j] = {
              comparisonMapping.createMeta(),
              .seedAStartPos = static_cast<int32_t>(cmp.seeds[j].seedAStartPos),
              .seedBStartPos = static_cast<int32_t>(cmp.seeds[j].seedBStartPos),
            };
            // PLOGE << "Seeds " << bucketMeta[i * NSEEDS + j].seedAStartPos << " " <<  bucketMeta[i * NSEEDS + j].seedBStartPos;
            // PLOGE << "Lengths " << bucketMeta[i * NSEEDS + j].sizeA << " " <<  bucketMeta[i * NSEEDS + j].sizeB;
            auto validCmp = bucketMeta[i * NSEEDS + j].seedAStartPos != -1;

            #pragma omp atomic
            cmps_total += validCmp;

            batch.origin_comparison_index[algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex + i * NSEEDS + j] = packOriginIndex(cmp.originalComparisonIndex, j);
            cellCount += cmp.sizeA * cmp.sizeB * validCmp;
          }
        } else {
          auto* bucketMeta = (SWMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
          bucketMeta[i] = comparisonMapping.createMeta();
          batch.origin_comparison_index[algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex + i] = packOriginIndex(cmp.originalComparisonIndex, 0);
          cellCount += cmp.sizeA * cmp.sizeB;
        }

        comparisonCount++;
      }
      cmpT.tock();
    }
  }
  stageTimers[2].tock();

  stageTimers[0].tock();
  json logData = {
    {"sequences_count", seqs.size() },
    {"comparisons_count", cmps_total },
    {"batches_created", batches.size()},
    {"time_total", stageTimers[0].seconds()},
    {"time_partition", stageTimers[1].seconds()},
    {"time_copy_inputs", stageTimers[2].seconds()},
    {"time_sequence", seqT.accumulate_microseconds() / 1e6},
    {"time_cmp", cmpT.accumulate_microseconds() / 1e6},
  };
  PLOGD << "BATCHCREATE: " << logData.dump();
  return batches;
}

template
std::vector<Batch> create_batches<Comparison>(const RawSequences& seqs, Comparisons& cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config);

template
std::vector<Batch> create_batches<MultiComparison>(const RawSequences& seqs, MultiComparisons& cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config);

}