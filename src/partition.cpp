#include "partition.h"

#include <tuple>

namespace ipu {
namespace partition {
  int fillFirst(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity) {
    mapping = std::vector<std::tuple<int, int>>(A.size(), {0, 0});
    std::vector<BucketData> buckets(bucketCount, {0, 0, 0, 0});
    int bucketIndex = 0;
    for (int i = 0; i < A.size(); ++i) {
      const auto& a = A[i];
      const auto& b = B[i];

      // find next empty bucket
      while (bucketIndex < bucketCount) {
        auto& [bN, bA, bB, _] = buckets[bucketIndex];
        if (bN + 1 <= bucketCountCapacity && bA + a.size() <= bucketCapacity && bB + b.size() <= bucketCapacity) {
          mapping[i] = {bucketIndex, i};
          bN++;
          bA += a.size();
          bB += b.size();
          // PLOGD << "i: " << i << " bucket index: " << bucketIndex << " cap: " << bA << "/" << bucketCapacity << " n: " << bN << "/" << bucketCountCapacity << "\n";
          break;
        } else {
          bucketIndex++;
        }
      }
      if (bucketIndex >= bucketCount) {
          return 1;
      }
    }

    return 0;
  }

  int roundRobin(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity) {
    mapping = std::vector<std::tuple<int, int>>(A.size(), {0, 0});
    std::vector<BucketData> buckets(bucketCount, {0, 0, 0, 0});
    int bucketIndex = 0;
    for (int i = 0; i < A.size(); ++i) {
      const auto& a = A[i];
      const auto& b = B[i];

      // find next empty bucket
      int boff = 0;
      for (; boff < bucketCount; ++boff) {
        int bi = (bucketIndex + boff) % bucketCount;
        auto& [bN, bA, bB, _] = buckets[bi];
        if (bN + 1 > bucketCountCapacity || bA + a.size() > bucketCapacity || bB + b.size() > bucketCapacity) {
          continue;
        } else {
          mapping[i] = {bucketIndex, i};
          bN++;
          bA += a.size();
          bB += b.size();
          break;
        }
      }
      if (boff >= bucketCount) {
          return 1;
      }
      bucketIndex = (bucketIndex + 1) % bucketCount;
    }
    return 0;
  }

  // Greedy approach in which we always put current sequence into one with lowest weight
  int greedy(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity) {
    mapping = std::vector<std::tuple<int, int>>(A.size(), {0, 0});
    std::vector<BucketData> buckets(bucketCount, {0, 0, 0, 0});
    for (int i = 0; i < A.size(); ++i) {
      const auto& a = A[i];
      const auto& b = B[i];

      auto weight = a.size() * b.size();
      int smallestBucket = -1;
      int smallestBucketWeight = 0;
      for (int bi = 0; bi < bucketCount; ++bi) {
        auto [bN, bA, bB, bW] = buckets[bi];
        if (!(bN + 1 > bucketCountCapacity || bA + a.size() > bucketCapacity || bB + b.size() > bucketCapacity)) {
          if (smallestBucket == -1 || smallestBucketWeight > bW) {
            smallestBucket = bi;
            smallestBucketWeight = bW;
          }
        }
      }

      if (smallestBucket == -1) {
        return 1;
      }

      auto& [bN, bA, bB, bW] = buckets[smallestBucket];
      bN++;
      bA += a.size();
      bB += b.size();
      bW += weight;
      mapping[i] = {smallestBucket, i};
    }
    return 0;
    // return mapping;
  }
}}