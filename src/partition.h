#ifndef PARTITION_H
#define PARTITION_H
#include <string>
#include <vector>
#include <queue>
#include <list>
#include "types.h"
#include "ipu_config.h"

namespace ipu { namespace partition {
  // mapping of the comparison to an individual tile
  struct ComparisonMapping {
    int comparisonIndex;  // index in cmp input list
    size_t sizeA;
    size_t offsetA;
    size_t sizeB;
    size_t offsetB;
    size_t seedAStartPos;
    size_t seedBStartPos;

    std::string toString() const;
  };

  struct SequenceMapping {
    int index;  // index in input list
    size_t offset;

    std::string toString() const;
  };

  struct ComparisonData {
    int comparisonIndex;
    int indexA;
    int indexB;
    size_t sizeA;
    size_t sizeB;
    int seedAStartPos;
    int seedBStartPos;

    size_t complexity;

    bool operator<(const ComparisonData& other) const;
  };

  struct Bucket {
    int bucketIndex;

    size_t sequenceCapacity;
    size_t comparisonCapacity;

    size_t longestLength;
    size_t totalSequenceLength;

    uint64_t totalCells;

    std::vector<SequenceMapping> seqs;
    std::vector<ComparisonMapping> cmps;

    Bucket(int bucketIndex, size_t sequenceCapacity, size_t comparisonCapacity);

    bool addComparison(int comparisonIndex, int indexA, int indexB, size_t sizeA, size_t sizeB, size_t seedAStartPos, size_t seedBStartPos);
    bool addComparison(const ComparisonData&);

    std::string toString() const;
  };

  bool operator<(const ipu::partition::Bucket& b1, const ipu::partition::Bucket& b2);
  bool operator>(const ipu::partition::Bucket& b1, const ipu::partition::Bucket& b2);

  struct BatchMapping {
    std::vector<Bucket> buckets;

    BatchMapping(int bucketCount, size_t sequenceCapacity, size_t comparisonCapacity);
    BatchMapping(IPUAlgoConfig config);

    std::string toString() const;
  };

  using BucketHeapRef = std::priority_queue<std::reference_wrapper<Bucket>, std::deque<std::reference_wrapper<Bucket>>, std::greater<std::deque<std::reference_wrapper<Bucket>>::value_type>>;
  using BucketHeap = std::priority_queue<Bucket, std::deque<Bucket>, std::greater<std::deque<Bucket>::value_type>>;

  // generic methods
  std::list<BatchMapping> mapBatches(IPUAlgoConfig config, const RawSequences& Seqs, const Comparisons& Cmps);
}}

#endif