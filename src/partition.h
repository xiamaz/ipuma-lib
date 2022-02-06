#ifndef PARTITION_H
#define PARTITION_H
#include<string>
#include<vector>

namespace ipu { namespace partition {
  enum class Algorithm {fillFirst, roundRobin, greedy};

  struct BucketData {
    int count;
    int lenA;
    int lenB;
    int weight;
  };

  int fillFirst(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity);
  int roundRobin(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity);
  int greedy(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity);
}}

#endif