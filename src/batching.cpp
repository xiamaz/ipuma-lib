#include "batching.h"
#include <plog/Log.h>
#include <omp.h>

#include "partition.h"

namespace ipu {
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

void Batch::initialize(IPUAlgoConfig config) {
  inputs.resize(config.getInputBufferSize32b(), 0);
  results.resize(config.getTotalNumberOfComparisons() * 3, 0);
  origin_comparison_index.resize(config.getTotalNumberOfComparisons(), -1);
  metaOffset = config.getOffsetMetadata();
  maxComparisons = config.getTotalNumberOfComparisons();
  numComparisons = 0;
  cellCount = 0;
  dataCount = 0;
}

BlockAlignmentResults Batch::get_result() {
  std::vector<std::array<int32_t, NSEEDS>> scores(maxComparisons);
  std::vector<std::array<uint32_t, NSEEDS>> a_range_result(maxComparisons);
  std::vector<std::array<uint32_t, NSEEDS>> b_range_result(maxComparisons);

  // TODO this should not be hardcoded here but passed as a config to make layout changes easier.
  int aOffset = results.size() / 3;
  int bOffset = aOffset * 2;

  // mapping is now responsibility of the caller
  for (int i = 0; i < maxComparisons; ++i) {
    int aindex = i + aOffset;
    int bindex = i + bOffset;
    for (size_t j = 0; j < NSEEDS; j++) {
      // if (results[i*NSEEDS + j] != 0) {
        // PLOGW << results[i*NSEEDS + j];
        // PLOGW <<"AAAAAAAAAAAAAAAAAAAAAAAA" << i << "   "<< j;
        // PLOGW.printf("scores[%d/NSEEDS][%d] = results[%d*NSEEDS + %d]\n", i, j, i, j);
      // }
      scores[i][j] = results[i*NSEEDS + j];
      a_range_result[i][j] = results[aOffset + i*NSEEDS + j];
      b_range_result[i][j] = results[bOffset + i*NSEEDS + j];
    }
    // scores[i] = results[i];
    // a_range_result[i] = results[aindex];
    // b_range_result[i] = results[bindex];
  }

  return {scores, a_range_result, b_range_result};
}

std::vector<Batch> create_batches(const RawSequences& seqs, Comparisons& cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config) {
  std::vector<swatlib::TickTock> stageTimers(3);
  stageTimers[0].tick();
  stageTimers[1].tick();
  auto mappings = partition::mapBatches(algoconfig, seqs, cmps);
  stageTimers[1].tock();

  std::vector<Batch> batches(mappings.size());
  const auto inputBufferSize = algoconfig.getInputBufferSize32b();
  auto encodeTable = swatlib::getEncoder(config.datatype).getCodeTable();
  const bool isSeeded = algoconfig.vtype == ipu::VertexType::xdropseedextend || algoconfig.vtype == VertexType::xdroprestrictedseedextend;

  stageTimers[2].tick();
  size_t vertexMetaSize = algoconfig.maxComparisonsPerVertex * sizeof(SWMeta);
  if (isSeeded) {
    vertexMetaSize = algoconfig.maxComparisonsPerVertex * sizeof(XDropMeta);
  }

  swatlib::TickTock seqT;
  swatlib::TickTock cmpT;
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
        auto* cmp = comparisonMapping.comparison;
        if (isSeeded) {
          auto* bucketMeta = (XDropMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
          for (int j = 0; j < NSEEDS; ++j) {
            bucketMeta[i * NSEEDS + j] = {
              comparisonMapping.createMeta(),
              .seedAStartPos = static_cast<int32_t>(cmp->seeds[j].seedAStartPos),
              .seedBStartPos = static_cast<int32_t>(cmp->seeds[j].seedBStartPos),
            };
          }
        } else {
          auto* bucketMeta = (SWMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
          bucketMeta[i] = comparisonMapping.createMeta();
        }

        batch.origin_comparison_index[algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex + i] = cmp->originalComparisonIndex;

        cellCount += cmp->sizeA * cmp->sizeB;
        comparisonCount++;
      }
      cmpT.tock();
    }
  }
  stageTimers[2].tock();

  auto cmps_total = 0;
  for (auto &&cmp : cmps) {
    for (auto &&seed : cmp.seeds) {
      cmps_total += seed.seedAStartPos != -1;
    }
  }

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

}