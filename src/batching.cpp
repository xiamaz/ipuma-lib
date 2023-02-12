#include "batching.h"
#include <plog/Log.h>

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

Batch::Batch(IPUAlgoConfig config) : inputs(config.getInputBufferSize32b(), 0), results(config.getTotalNumberOfComparisons() * 3, 0), origin_comparison_index(config.getTotalNumberOfComparisons(), -1), metaOffset(config.getOffsetMetadata()), numComparisons{0}, cellCount{0}, dataCount{0} {
}

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

std::vector<Batch> create_batches(const RawSequences& seqs, const Comparisons& cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config) {
  std::vector<swatlib::TickTock> stageTimers(3);
  stageTimers[0].tick();
  stageTimers[1].tick();
  auto mappings = partition::mapBatches(algoconfig, seqs, cmps);
  stageTimers[1].tock();

  std::vector<Batch> batches;
  const auto inputBufferSize = algoconfig.getInputBufferSize32b();
  auto encodeTable = swatlib::getEncoder(config.datatype).getCodeTable();
  const bool isSeeded = algoconfig.vtype == ipu::VertexType::xdropseedextend;
  

  stageTimers[2].tick();
  size_t vertexMetaSize = algoconfig.maxComparisonsPerVertex * sizeof(SWMeta);
  if (isSeeded) {
    vertexMetaSize = algoconfig.maxComparisonsPerVertex * sizeof(XDropMeta);
  }
  for (const auto& map : mappings) {
    batches.push_back({algoconfig});
    Batch& batch = batches.back();

    int8_t* seqInput = batch.getSequenceBuffer();
    int8_t* metaInput = batch.getMetaBuffer();
    auto& cellCount = batch.cellCount;
    auto& dataCount = batch.dataCount;
    auto& comparisonCount = batch.numComparisons;

    for (const auto& bucketMapping : map.buckets) {
      const size_t offsetSequence = bucketMapping.bucketIndex * algoconfig.getBufsize32b() * 4;

      auto* bucketSeq = seqInput + offsetSequence;
      for (const auto& sequenceMapping : bucketMapping.seqs) {
        const char *seq = seqs[sequenceMapping.index].data();
        size_t seqSize = seqs[sequenceMapping.index].size();
        for (int j = 0; j < seqSize; ++j) {
          bucketSeq[sequenceMapping.offset + j] = encodeTable[seq[j]];
        }
        dataCount += seqSize;
      }

      for (int i = 0; i < bucketMapping.cmps.size(); ++i) {
        const auto& comparison = bucketMapping.cmps[i];
        if (isSeeded) {
          auto* bucketMeta = (XDropMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
          bucketMeta[i] = {
            {
              .sizeA = static_cast<int32_t>(comparison.sizeA),
              .offsetA = static_cast<int32_t>(comparison.offsetA),
              .sizeB = static_cast<int32_t>(comparison.sizeB),
              .offsetB = static_cast<int32_t>(comparison.offsetB),
            },
            .seedAStartPos = static_cast<int32_t>(comparison.seedAStartPos),
            .seedBStartPos = static_cast<int32_t>(comparison.seedBStartPos),
          };
        } else {
          auto* bucketMeta = (SWMeta*)(metaInput) + algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex;
          bucketMeta[i] = {
            .sizeA = static_cast<int32_t>(comparison.sizeA),
            .offsetA = static_cast<int32_t>(comparison.offsetA),
            .sizeB = static_cast<int32_t>(comparison.sizeB),
            .offsetB = static_cast<int32_t>(comparison.offsetB),
          };
        }

        batch.origin_comparison_index[algoconfig.maxComparisonsPerVertex * bucketMapping.bucketIndex + i] = comparison.comparisonIndex;

        cellCount += comparison.sizeA * comparison.sizeB;
        comparisonCount++;
      }
    }
  }
  stageTimers[2].tock();


  stageTimers[0].tock();
  json logData = {
    {"sequences_count", seqs.size() },
    {"comparisons_count", cmps.size() },
    {"batches_created", batches.size()},
    {"time_total", stageTimers[0].seconds()},
    {"time_partition", stageTimers[1].seconds()},
    {"time_copy_inputs", stageTimers[2].seconds()},
  };
  PLOGD << "BATCHCREATE: " << logData.dump();
  return batches;
}

}