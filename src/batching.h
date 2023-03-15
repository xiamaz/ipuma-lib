#pragma once

#include <vector>
#include "ipu_base.h"
#include "types.h"

namespace ipu {

typedef int32_t OriginIndex;

struct BlockAlignmentResults {
  std::vector<int32_t> scores;
  std::vector<uint32_t> a_range_result;
  std::vector<uint32_t> b_range_result;
};

int32_t getComparisonIndex(OriginIndex);
int32_t getSeedIndex(OriginIndex);
std::tuple<int32_t, int32_t> unpackOriginIndex(OriginIndex);

struct Batch {
  std::vector<int32_t> inputs;
  std::vector<int32_t> results;
  std::vector<OriginIndex> origin_comparison_index;

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

template<typename C>
std::vector<Batch> create_batches(const RawSequences& seqs, std::vector<C>& cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config);

extern template
std::vector<Batch> create_batches<Comparison>(const RawSequences& Seqs, Comparisons& Cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config);

extern template
std::vector<Batch> create_batches<MultiComparison>(const RawSequences& Seqs, MultiComparisons& Cmps, const IPUAlgoConfig& algoconfig, const SWConfig& config);
}