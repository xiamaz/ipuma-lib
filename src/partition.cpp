#include "partition.h"

#include <sstream>
#include <functional>
#include <queue>
#include <tuple>
#include <algorithm>
#include <array>
#include <iostream>
#include <math.h>
#include <limits.h>

#include <plog/Log.h>


#include <iostream>
#include <omp.h>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace ipu {
namespace partition {
  std::string ComparisonMapping::toString() const {
    std::stringstream ss;
    ss << "CM[" << comparisonIndex << ": a(l" << sizeA << " o" << offsetA << ") b(l" << sizeB << " o" << offsetB << ")]";
    return ss.str();
  }

  std::string SequenceMapping::toString() const {
    std::stringstream ss;
    ss << "SM[" << index << ": o" << offset << "]";
    return ss.str();
  }

  Bucket::Bucket(int bucketIndex, size_t sequenceCapacity, size_t comparisonCapacity) : bucketIndex(bucketIndex), sequenceCapacity(sequenceCapacity), comparisonCapacity(comparisonCapacity) {
    longestLength = 0;
    totalSequenceLength = 0;
  }

  std::string Bucket::toString() const {
    std::stringstream ss;
    ss << "BMapping[" << bucketIndex << ": cmps(" << cmps.size() << ") maxLen(" << longestLength << ") seqSize(" << totalSequenceLength << ")]";
    return ss.str();
  }

  bool Bucket::addComparison(int comparisonIndex, int indexA, int indexB, size_t sizeA, size_t sizeB) {
    size_t newTotalLength = totalSequenceLength + sizeA + sizeB;
    size_t newTotalCmps = cmps.size() + 1;

    // PLOGW << "TL: " << newTotalLength;
    // PLOGW << sequenceCapacity;
    // PLOGW << newTotalCmps;
    // PLOGW << comparisonCapacity;

    if (newTotalLength < sequenceCapacity && newTotalCmps < comparisonCapacity) {
      size_t offsetA = totalSequenceLength;
      size_t offsetB = totalSequenceLength + sizeA;

      seqs.push_back({
        .index = indexA,
        .offset = offsetA,
      });
      seqs.push_back({
        .index = indexB,
        .offset = offsetB,
      });

      cmps.push_back({
        .comparisonIndex = comparisonIndex,
        .sizeA = sizeA,
        .offsetA = offsetA,
        .sizeB = sizeB,
        .offsetB = offsetB,
      });

      totalSequenceLength = newTotalLength;
      longestLength = std::max({longestLength, sizeA, sizeB});
      return true;
    }
    return false;
  }

  std::string BatchMapping::toString() const {
    std::stringstream ss;
    ss << "BMap[" << buckets.size() << "]";
    return ss.str();
  }

  BatchMapping::BatchMapping(int bucketCount, size_t sequenceCapacity, size_t comparisonCapacity) {
    for (int i = 0; i < bucketCount; ++i) {
      buckets.push_back(Bucket(i, sequenceCapacity, comparisonCapacity));
    }
  }

  BatchMapping::BatchMapping(IPUAlgoConfig config) {
    for (int i = 0; i < config.tilesUsed; ++i) {
      buckets.push_back(Bucket(i, config.bufsize, config.maxBatches));
    }
  }

  struct ComparisonData {
    int comparisonIndex;
    int indexA;
    int indexB;
    size_t sizeA;
    size_t sizeB;
  };

  bool batchFillFirst(BatchMapping& map, const ComparisonData& cmpData) {
    for (auto& bucket : map.buckets) {
      if (bucket.addComparison(cmpData.comparisonIndex, cmpData.indexA, cmpData.indexB, cmpData.sizeA, cmpData.sizeB)) {
        return true;
      }
    }
    return false;
  }

  bool batchRoundRobin(BatchMapping& map, const ComparisonData& cmpData) {
    return false;
  }
  bool batchGreedy(BatchMapping& map, const ComparisonData& cmpData) {
    return false;
  }

  using MappingFunction = std::function<bool(BatchMapping&, const ComparisonData&)>;

  MappingFunction selectAlgorithm(Algorithm algo) {
    switch (algo) {
    case Algorithm::fillFirst:
      return batchFillFirst;
      break;
    case Algorithm::roundRobin:
      return batchRoundRobin;
      break;
    case Algorithm::greedy:
      return batchGreedy;
      break;
    }
    throw std::runtime_error("No valid algorithm given");
  }

  ComparisonData createCmpData(int i, const Comparison& cmp, const RawSequences& seqs) {
    return {
      .comparisonIndex = i,
      .indexA = cmp.indexA,
      .indexB = cmp.indexB,
      .sizeA = seqs[cmp.indexA].size(),
      .sizeB = seqs[cmp.indexB].size(),
    };
  }

  std::list<BatchMapping> mapBatches(IPUAlgoConfig config, const RawSequences& Seqs, const Comparisons& Cmps) {
    std::list<BatchMapping> mappings;

    MappingFunction mapData = selectAlgorithm(config.fillAlgo);

    mappings.push_back(BatchMapping(config));
    BatchMapping& curBatch = mappings.back();
    for (int i = 0; i < Cmps.size(); ++i) {
      const Comparison& cmp = Cmps[i];
      ComparisonData cmpData = createCmpData(i, cmp, Seqs);
      if (!mapData(curBatch, cmpData)) {
        mappings.push_back(BatchMapping(config));
        curBatch = mappings.back();
        if (!mapData(curBatch, cmpData)) {
          throw std::runtime_error("Could not insert data into a new batch");
        }
      }
    }
    return mappings;
  }
}}