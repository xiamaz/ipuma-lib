#ifndef PARTITION_H
#define PARTITION_H
#include<string>
#include<vector>
#include<queue>
#include "types.h"

namespace ipu { namespace partition {
  enum class SequenceOrigin{unordered, A, B};
  std::string sequenceOriginToString(SequenceOrigin o);

  struct ComparisonMapping {
    int comparisonIndex;  // index in cmp input list
    int sizeA;
    int offsetA;
    int sizeB;
    int offsetB;

    std::string toString();
  };

  struct SequenceMapping {
    int index;  // index in input list
    int offset;
    SequenceOrigin origin = SequenceOrigin::unordered;

    std::string toString();
  };

  struct BucketMapping {
    int bucketIndex;
    int maxLen; // length of longest sequence
    int seqSize; // length of current sequences
    int weight; // custom weight (currently only used in round robin)
    std::vector<SequenceMapping> seqs;
    std::vector<ComparisonMapping> cmps;

    std::string toString();
  };

  bool operator<(const ipu::partition::BucketMapping& b1, const ipu::partition::BucketMapping& b2);
  bool operator>(const ipu::partition::BucketMapping& b1, const ipu::partition::BucketMapping& b2);

  struct BucketMap {
    std::vector<BucketMapping> buckets;
    int numBuckets;
    int cmpCapacity;
    int sequenceCapacity;

    BucketMap();
    BucketMap(int nB, int nC, int sC);

    std::string toString();
  };

  struct BucketData {
    int count;
    int lenSeq;
    int weight;
  };

  using BucketHeap = std::priority_queue<std::reference_wrapper<BucketMapping>, std::deque<std::reference_wrapper<BucketMapping>>, std::greater<std::deque<std::reference_wrapper<BucketMapping>>::value_type>>;

  bool fillFirst(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, int& curBucket);
  bool roundRobin(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, int& curBucket);
  bool greedy(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, BucketHeap& heap);

  void fillFirst(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset = 0);
  void roundRobin(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset = 0);
  void greedy(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset = 0);

  void fillFirst(BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int indexOffset = 0);
  void roundRobin(BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int indexOffset = 0);
  void greedy(BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int indexOffset = 0);

  // generic methods
  bool fillBuckets(Algorithm algo, BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, int& curBucket, BucketHeap& heap);
}}

#endif