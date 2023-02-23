#pragma once

#include <vector>
#include "ipu_base.h"

namespace ipu {
struct BlockAlignmentResults {
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
};

struct Batch {
  std::vector<int32_t> inputs;
  std::vector<int32_t> results;
  std::vector<int32_t> origin_comparison_index;

  Batch();
  Batch(IPUAlgoConfig config);

  void initialize(IPUAlgoConfig config);

  size_t maxComparisons;
  size_t numComparisons;
  size_t metaOffset;

  // metrics kept for computing performance
  uint64_t cellCount;
  uint64_t dataCount;

  swatlib::TickTock tick;

  std::string toString() const;

  BlockAlignmentResults get_result();

  int8_t* getSequenceBuffer();
  int8_t* getMetaBuffer();
};

std::vector<Batch> create_batches(const RawSequences& Seqs, const Comparisons& Cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config);
}