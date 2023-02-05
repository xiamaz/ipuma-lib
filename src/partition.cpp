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
      totalCells += sizeA * sizeB;
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

  class PartitioningAlgorithm {
  protected:
    const IPUAlgoConfig& config;

    virtual BatchMapping& createBatch() = 0;
  public:
    std::list<BatchMapping> mappings;

    PartitioningAlgorithm(const IPUAlgoConfig& config) : config(config), mappings{BatchMapping(config)} {
    }
    virtual bool addComparison(const ComparisonData& cmpData, BatchMapping& curMapping) = 0;

    std::list<BatchMapping> getMappings() {
      return std::move(mappings);
    }

    bool addComparison(const ComparisonData& cmpData) {
      BatchMapping& curMapping = mappings.back();
      if (!addComparison(cmpData, curMapping)) {
        curMapping = createBatch();
        if (!addComparison(cmpData, curMapping)) {
          throw std::runtime_error("Could not insert data into a new batch");
        }
      }
      return true;
    }
  };

  class FillFirstPartitioning : public PartitioningAlgorithm {
    using PartitioningAlgorithm::PartitioningAlgorithm;

    int bucketIndex = 0;

    BatchMapping& createBatch() override {
      mappings.push_back(BatchMapping(config));
      bucketIndex = 0;
      return mappings.back();
    }

    bool addComparison(const ComparisonData& cmpData, BatchMapping& curMapping) override {
      for (; bucketIndex < curMapping.buckets.size(); ++bucketIndex) {
        Bucket& bucket = curMapping.buckets[bucketIndex];
        if (bucket.addComparison(cmpData.comparisonIndex, cmpData.indexA, cmpData.indexB, cmpData.sizeA, cmpData.sizeB)) {
          return true;
        }
      }
      return false;
    }
  };

  using BucketWrapper = std::reference_wrapper<Bucket>;

  class BucketCompare {
  public:
    bool operator() (const Bucket& a, const Bucket& b) {
      // auto a_remaining = a.sequenceCapacity - a.totalSequenceLength;
      // auto b_remaining = b.sequenceCapacity - b.totalSequenceLength;
      // return a_remaining <= b_remaining;
      return a.totalCells > b.totalCells;
    }
  };

  class GreedyPartitioning : public PartitioningAlgorithm {
    virtual bool addComparison(const ComparisonData& cmpData, BatchMapping& curMapping) override {
      BatchMapping& m = mappings.back();
      std::pop_heap(m.buckets.begin(), m.buckets.end(), cmp);
      Bucket& b = m.buckets.back();
      auto success = b.addComparison(cmpData.comparisonIndex, cmpData.indexA, cmpData.indexB, cmpData.sizeA, cmpData.sizeB);
      std::push_heap(m.buckets.begin(), m.buckets.end(), cmp);
      return success;
    }

    BatchMapping& createBatch() override {
      mappings.push_back(BatchMapping(config));
      BatchMapping& m = mappings.back();
      std::make_heap(m.buckets.begin(), m.buckets.end(), cmp);
      return mappings.back();
    }

  public:
    BucketCompare cmp;
    GreedyPartitioning(const IPUAlgoConfig& config) : PartitioningAlgorithm(config), cmp{} {
      BatchMapping& m = mappings.back();
      std::make_heap(m.buckets.begin(), m.buckets.end(), BucketCompare());
    }
  };

  std::unique_ptr<PartitioningAlgorithm> createPartitioner(const IPUAlgoConfig& config) {
    switch (config.fillAlgo) {
    case Algorithm::fillFirst:
      return std::unique_ptr<PartitioningAlgorithm>(new FillFirstPartitioning(config));
      break;
    case Algorithm::roundRobin:
      throw std::runtime_error("Not implemented");
      break;
    case Algorithm::greedy:
      return std::unique_ptr<PartitioningAlgorithm>(new GreedyPartitioning(config));
      break;
    }
    throw std::runtime_error("Not any of supported algorithm");
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

    PLOGD << "Mapping " << Cmps.size() << " comparisons to batches";

    auto Partitioner = createPartitioner(config);
    for (int i = 0; i < Cmps.size(); ++i) {
      ComparisonData cmpData = createCmpData(i, Cmps[i], Seqs);
      Partitioner->addComparison(cmpData);
    }
    PLOGD << "Created " << Partitioner->mappings.size() << " batches";
    return Partitioner->getMappings();
  }
}}