#include "partition.h"

#include <sstream>
#include <functional>
#include <queue>
#include <tuple>
#include <algorithm>
#include <array>
#include <iostream>
#include <math.h>

#include <plog/Log.h>


#include <iostream>

namespace ipu {
namespace partition {
  std::string sequenceOriginToString(SequenceOrigin o) {
    switch(o) {
    case SequenceOrigin::A:
      return "SequenceOrigin::A";
    case SequenceOrigin::B:
      return "SequenceOrigin::B";
    case SequenceOrigin::unordered:
      return "SequenceOrigin::unordered";
    }
    throw std::runtime_error("Unhandled sequence type.");
  }

  std::string ComparisonMapping::toString() {
    std::stringstream ss;
    ss << "CM[" << comparisonIndex << ": a(l" << sizeA << " o" << offsetA << ") b(l" << sizeB << " o" << offsetB << ")]";
    return ss.str();
  }

  std::string SequenceMapping::toString() {
    std::stringstream ss;
    ss << "SM[" << index << ": o" << offset << " t" << sequenceOriginToString(origin) << "]";
    return ss.str();
  }

  std::string BucketMapping::toString() {
    std::stringstream ss;
    ss << "BMapping[" << bucketIndex << ": cmps(" << cmps.size() << ") maxLen(" << maxLen << ") seqSize(" << seqSize << ") weight(" << weight << ")]";
    return ss.str();
  }

  std::string BucketMap::toString() {
    std::stringstream ss;
    ss << "BMap[" << numBuckets << ": maxCmps(" << cmpCapacity << ") bufsize(" << sequenceCapacity << ")]";
    return ss.str();
  }

  BucketMap::BucketMap() : numBuckets(0), cmpCapacity(0), sequenceCapacity(0) { }

  BucketMap::BucketMap(int nB, int cC, int sC) : numBuckets(nB), cmpCapacity(cC), sequenceCapacity(sC) {
    buckets.resize(nB);
    for (int i = 0; i < nB; ++i) {
      buckets[i].bucketIndex = i;
    }
  }

  bool operator>(const ipu::partition::BucketMapping& b1, const ipu::partition::BucketMapping& b2) {
    return b1.weight > b2.weight;
  }

  bool operator<(const ipu::partition::BucketMapping& b1, const ipu::partition::BucketMapping& b2) {
    return b1.weight < b2.weight;
  }

  struct SequenceInfo {
    const int aLen;
    const int bLen;
    int effALen;
    int effBLen;
    int offsetA;
    int offsetB;
  };

  BucketMapping& getNextBucket(BucketMap& map, int curBucket, const Comparison& cmp, SequenceInfo& info, int indexOffset) {
      int offset = 0;
      int effAlen, offsetA, effBlen, offsetB;
      for (; offset < map.numBuckets; ++offset) {
        int wrapped = (curBucket + offset) % map.numBuckets;
        const auto& bucket = map.buckets[wrapped];

        const SequenceMapping* smap = nullptr;
        effAlen = info.aLen;
        effBlen = info.bLen;
        for (int n = 0; n < bucket.seqs.size(); ++n) {
          const auto& s = bucket.seqs[n];
          if (s.index == indexOffset + cmp.indexA) {
            effAlen = 0;
            offsetA = s.offset;
          }
          if (s.index == indexOffset + cmp.indexB) {
            effBlen = 0;
            offsetB = s.offset;
          }
        }
        if ((bucket.cmps.size() + 1 <= map.cmpCapacity) && (bucket.seqSize + effAlen + effBlen <= map.sequenceCapacity)) {
          break;
        }
      }
      if (offset >= map.numBuckets) {
        throw std::runtime_error("Out of buckets.");
      }
      curBucket = (curBucket + offset) % map.numBuckets;
      info.effALen = effAlen;
      info.effBLen = effBlen;
      info.offsetA = offsetA;
      info.offsetB = offsetB;
      return map.buckets[curBucket];
  }

  void addCmpToBucket(BucketMapping& bucket, SequenceInfo& info, const Comparison& cmp, int bucketCmpIndex, int indexOffset) {
      bucket.maxLen = std::max(bucket.maxLen, static_cast<int>(info.aLen > info.bLen ? info.aLen : info.bLen));
      auto ai = indexOffset + cmp.indexA;
      auto bi = indexOffset + cmp.indexB;
      auto ci = indexOffset + bucketCmpIndex;
      if (info.effALen) {
        info.offsetA = bucket.seqSize;
        bucket.seqs.push_back({.index=ai, .offset=info.offsetA, .origin = SequenceOrigin::unordered});
      }
      if (info.effBLen) {
        info.offsetB = bucket.seqSize + info.effALen;
        bucket.seqs.push_back({.index=bi, .offset=info.offsetB, .origin = SequenceOrigin::unordered});
      }
      bucket.cmps.push_back({.comparisonIndex = ci, .sizeA = static_cast<int>(info.aLen), .offsetA = info.offsetA, .sizeB = static_cast<int>(info.bLen), .offsetB = info.offsetB});
      bucket.seqSize += info.effALen + info.effBLen;
  }

  void fillFirst(BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int indexOffset) {
    int curBucket = 0;
    for (int i = 0; i < Cmps.size(); ++i) {
      const auto& cmp = Cmps[i];
      SequenceInfo info = {
        .aLen = static_cast<int>(Seqs[cmp.indexA].size()),
        .bLen = static_cast<int>(Seqs[cmp.indexB].size()),
        0, 0, 0, 0
      };

      auto& bucket = getNextBucket(map, curBucket, cmp, info, indexOffset);
      curBucket = bucket.bucketIndex;
      addCmpToBucket(bucket, info, cmp, i, indexOffset);
    }
  }

  void fillFirst(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset) {
    int curBucket = 0;
    if (!fillFirst(map, A, B, indexOffset, curBucket)) {
      throw std::runtime_error("Out of buckets.");
    }
  }
  void roundRobin(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset) {
    int curBucket = 0;
    if (!roundRobin(map, A, B, indexOffset, curBucket)){
      throw std::runtime_error("Out of buckets.");
    }
  }
  void greedy(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset) {
    BucketHeap q;
    for (auto& b : map.buckets) {
      q.push(std::ref(b));
    }
    if (!greedy(map, A, B, indexOffset, q)) {
      throw std::runtime_error("Out of buckets.");
    }
  }


  bool fillFirst(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, int& curBucket) {
    for (int i = 0; i < A.size(); ++i) {
      const auto aLen = A[i].size();
      const auto bLen = B[i].size();

      int offset = 0;
      for (; offset < map.numBuckets; ++offset) {
        int wrapped = (curBucket + offset) % map.numBuckets;
        const auto& bucket = map.buckets[wrapped];
        if ((bucket.cmps.size() + 1 <= map.cmpCapacity) && (bucket.seqSize + aLen + bLen <= map.sequenceCapacity)) break;
      }
      if (offset >= map.numBuckets) {
        return false;
      }
      curBucket = (curBucket + offset) % map.numBuckets;
      auto& bucket = map.buckets[curBucket];
      bucket.maxLen = std::max(bucket.maxLen, static_cast<int>(aLen > bLen ? aLen : bLen));
      auto ci = indexOffset + i;
      int offsetA = bucket.seqSize;
      int offsetB = bucket.seqSize + aLen;
      bucket.seqs.push_back({.index=ci, .offset=offsetA, .origin = SequenceOrigin::A});
      bucket.seqs.push_back({.index=ci, .offset=offsetB, .origin = SequenceOrigin::B});
      bucket.cmps.push_back({.comparisonIndex = ci, .sizeA = static_cast<int>(aLen), .offsetA = offsetA, .sizeB = static_cast<int>(bLen), .offsetB = offsetB});
      bucket.seqSize += aLen + bLen;
    }
    return true;
  }

  void roundRobin(BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int indexOffset) {
    int curBucket = 0;
    for (int i = 0; i < Cmps.size(); ++i) {
      const auto& cmp = Cmps[i];
      SequenceInfo info = {
        .aLen = static_cast<int>(Seqs[cmp.indexA].size()),
        .bLen = static_cast<int>(Seqs[cmp.indexB].size()),
        0, 0, 0, 0
      };

      auto& bucket = getNextBucket(map, curBucket, cmp, info, indexOffset);
      addCmpToBucket(bucket, info, cmp, i, indexOffset);
      curBucket = bucket.bucketIndex + 1;
    }
  }

  bool roundRobin(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, int& curBucket) {
    for (int i = 0; i < A.size(); ++i) {
      const auto& aLen = A[i].size();
      const auto& bLen = B[i].size();

      int offset = 0;
      for (; offset < map.numBuckets; ++offset) {
        int wrapped = (curBucket + offset) % map.numBuckets;
        const auto& bucket = map.buckets[wrapped];
        if ((bucket.cmps.size() + 1 <= map.cmpCapacity) && (bucket.seqSize + aLen + bLen <= map.sequenceCapacity)) break;
      }
      if (offset >= map.numBuckets) {
        return false;
      }
      curBucket = (curBucket + offset) % map.numBuckets;
      auto& bucket = map.buckets[curBucket];
      bucket.maxLen = std::max(bucket.maxLen, static_cast<int>(aLen > bLen ? aLen : bLen));
      auto ci = indexOffset + i;
      int offsetA = bucket.seqSize;
      int offsetB = bucket.seqSize + aLen;
      bucket.seqs.push_back({.index=ci, .offset=offsetA, .origin = SequenceOrigin::A});
      bucket.seqs.push_back({.index=ci, .offset=offsetB, .origin = SequenceOrigin::B});
      bucket.cmps.push_back({.comparisonIndex = ci, .sizeA = static_cast<int>(aLen), .offsetA = offsetA, .sizeB = static_cast<int>(bLen), .offsetB = offsetB});
      bucket.seqSize += aLen + bLen;

      curBucket++; // increment current bucket for round robin
    }
    return true;
  }

  void greedy(BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int indexOffset) {
    // std::priority_queue<BucketMapping, std::deque<BucketMapping>, std::greater<BucketMapping&>> q;
    std::priority_queue<std::reference_wrapper<BucketMapping>, std::deque<std::reference_wrapper<BucketMapping>>, std::greater<std::deque<std::reference_wrapper<BucketMapping>>::value_type>> q;
    for (auto& b : map.buckets) {
      q.push(std::ref(b));
    }

    std::vector<std::pair<int, int>> srts(Cmps.size());
    for (size_t i = 0; i < Cmps.size(); i++) {
      const auto& cmp = Cmps[i];
      const auto aLen = Seqs[cmp.indexA].size();
      const auto bLen = Seqs[cmp.indexB].size();
      srts[i] = {aLen * bLen, i};
    }
    std::sort(srts.begin(), srts.end(), [](std::pair<int, int> &a, std::pair<int, int> &b) {
        return std::get<0>(a) > std::get<0>(b);
    });

    for (int i = 0; i < Cmps.size(); ++i) {
      const auto cmpIndex = std::get<1>(srts[i]);
      const auto& cmp = Cmps[cmpIndex];
      const auto aLen = Seqs[cmp.indexA].size();
      const auto bLen = Seqs[cmp.indexB].size();

      auto& bucket = q.top().get();
      q.pop();
      std::deque<std::reference_wrapper<BucketMapping>> qq;
      int tries = 0;
      int effALen, effBLen;
      int offsetA = 0;
      int offsetB = 0;
      while (tries < map.numBuckets) {
        effALen = aLen;
        effBLen = bLen;
        for (int n = 0; n < bucket.seqs.size(); ++n) {
          if (bucket.seqs[n].index == cmp.indexA) {
            effALen = 0;
            offsetA = bucket.seqs[n].offset;
          }
          if (bucket.seqs[n].index == cmp.indexB) {
            effBLen = 0;
            offsetB = bucket.seqs[n].offset;
          }
        }
        if ((bucket.cmps.size() + 1 > map.cmpCapacity) || (bucket.seqSize + effALen + effBLen > map.sequenceCapacity)) {
          qq.push_back(std::ref(bucket));
          bucket = q.top().get();
          q.pop();
          tries++;
          continue;
        } else {
          break;
        }
      }

      if (tries == map.numBuckets) {
          std::cout << bucket.toString() << "\n";
          for (auto& b : map.buckets) {
            std::cout << b.toString() << "\n";
          }
          throw std::runtime_error("Out of buckets");
      }

      for (auto b : qq) {
        q.push(b);
      }

      auto ai = indexOffset + cmp.indexA;
      auto bi = indexOffset + cmp.indexB;
      auto ci = indexOffset + cmpIndex;
      if (effALen) {
        offsetA = bucket.seqSize;
        bucket.seqs.push_back({.index=ai, .offset=offsetA, .origin = SequenceOrigin::unordered});
      }
      if (effBLen) {
        offsetB = bucket.seqSize + effALen;
        bucket.seqs.push_back({.index=bi, .offset=offsetB, .origin = SequenceOrigin::unordered});
      }
      bucket.cmps.push_back({.comparisonIndex = ci, .sizeA = static_cast<int>(aLen), .offsetA = offsetA, .sizeB = static_cast<int>(bLen), .offsetB = offsetB});

      bucket.maxLen = std::max(bucket.maxLen, static_cast<int>(aLen > bLen ? aLen : bLen));
      bucket.seqSize += effALen + effBLen;
      bucket.weight += aLen * bLen;
      q.push(std::ref(bucket));
    }
  }

  bool greedy(BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, BucketHeap& q) {
    std::vector<std::pair<int, int>> srts(A.size());
    for (size_t i = 0; i < A.size(); i++) {
      const auto aLen = A[i].size();
      const auto bLen = B[i].size();
      srts[i] = {aLen * bLen, i};
    }
    std::sort(srts.begin(), srts.end(), [](std::pair<int, int> &a, std::pair<int, int> &b) { 
      return std::get<0>(a) > std::get<0>(b);
    });

    for (int i = 0; i < A.size(); ++i) {
      const auto seqIndex = std::get<1>(srts[i]);
      const auto aLen = A[seqIndex].size();
      const auto bLen = B[seqIndex].size();

      if (aLen * bLen == 0) {
        PLOGE << "aLen " << aLen << " bLen " << bLen;
      }

      auto& bucket = q.top().get();
      q.pop();
      std::deque<std::reference_wrapper<BucketMapping>> qq;
      int tries = 0;
      while (tries < map.numBuckets) {
        if ((bucket.cmps.size() + 1 > map.cmpCapacity) || (bucket.seqSize + aLen + bLen > map.sequenceCapacity)) {
          qq.push_back(std::ref(bucket));
          bucket = q.top().get();
          q.pop();
          tries++;
          continue;
        } else {
          break;
        }
      }

      if (tries == map.numBuckets) {
        for (auto& b : map.buckets) {
          PLOGE << b.bucketIndex << ": " << map.sequenceCapacity - b.seqSize;
          PLOGE << "\t" << b.cmps.size() << "/" << map.cmpCapacity;
        }
        PLOGE << "Alen: " << aLen << " Blen: " << bLen;
        return false;
      }

      for (auto b : qq) {
        q.push(b);
      }

      auto ci = indexOffset + seqIndex;
      int offsetA = bucket.seqSize;
      int offsetB = bucket.seqSize + aLen;
      bucket.seqs.push_back({.index=ci, .offset=offsetA, .origin = SequenceOrigin::A});
      bucket.seqs.push_back({.index=ci, .offset=offsetB, .origin = SequenceOrigin::B});
      bucket.cmps.push_back({.comparisonIndex = ci, .sizeA = static_cast<int>(aLen), .offsetA = offsetA, .sizeB = static_cast<int>(bLen), .offsetB = offsetB});

      bucket.maxLen = std::max(bucket.maxLen, static_cast<int>(aLen > bLen ? aLen : bLen));
      bucket.seqSize += aLen + bLen;
      bucket.weight += aLen * bLen;
      q.push(std::ref(bucket));

    }
    return true;
  }

  bool fillBuckets(Algorithm algo, BucketMap& map, const RawSequences& A, const RawSequences& B, int indexOffset, int& curBucket, BucketHeap& heap) {
    bool ret = false;
    switch (algo) {
    case Algorithm::fillFirst:
      ret = partition::fillFirst(map, A, B, indexOffset, curBucket);
      break;
    case Algorithm::roundRobin:
      ret = partition::roundRobin(map, A, B, indexOffset, curBucket);
      break;
    case Algorithm::greedy:
      ret = partition::greedy(map, A, B, indexOffset, heap);
      break;
    }
    return ret;
  }
}}