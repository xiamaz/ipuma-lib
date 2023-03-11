#include "partition.h"
#include "complexity.h"

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
    ss << "CM[" << comparison->toString() << ": a(o" << offsetA << ") b(o" << offsetB << ")]";
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

  bool Bucket::addComparison(const Comparison& d) {
    size_t newTotalLength = totalSequenceLength + d.sizeA + d.sizeB;
    size_t newTotalCmps = cmps.size() * NSEEDS + NSEEDS;

    if (newTotalLength <= sequenceCapacity && newTotalCmps <= comparisonCapacity) {
      size_t offsetA = totalSequenceLength;
      size_t offsetB = totalSequenceLength + d.sizeA;

      seqs.push_back({
        .index = d.indexA,
        .offset = offsetA,
      });
      seqs.push_back({
        .index = d.indexB,
        .offset = offsetB,
      });

      if (d.seeds.size() != NSEEDS) {
        throw std::runtime_error("Number of seeds needs to match NSEEDS");
      }
      cmps.push_back({
        .offsetA = (size_t) offsetA,
        .offsetB = (size_t) offsetB,
        .comparison = &d
      });

      totalSequenceLength = newTotalLength;
      longestLength = std::max({longestLength, (size_t) d.sizeA, (size_t) d.sizeB});
      totalCells += d.sizeA * d.sizeB * NSEEDS;
      return true;
    }
    return false;
  }

  std::string BatchMapping::toString() const {
    std::stringstream ss;
    uint64_t totalSequenceSize = 0;
    uint64_t sequenceCapacity = 0;
    int totalCmps = 0;
    for (const auto& b : buckets) {
      totalCmps += b.cmps.size();
      sequenceCapacity += b.sequenceCapacity;
      totalSequenceSize += b.totalSequenceLength;
    }
    ss << "BMap[" << totalCmps << " cmps " << buckets.size() << " buckets] " << totalSequenceSize << "/" << sequenceCapacity;
    return ss.str();
  }

  BatchMapping::BatchMapping(int bucketCount, size_t sequenceCapacity, size_t comparisonCapacity) {
    for (int i = 0; i < bucketCount; ++i) {
      buckets.push_back(Bucket(i, sequenceCapacity, comparisonCapacity));
    }
  }

  BatchMapping::BatchMapping(IPUAlgoConfig config) {
    for (int i = 0; i < config.numVertices; ++i) {
      buckets.push_back(Bucket(i, config.vertexBufferSize, config.maxComparisonsPerVertex));
    }
  }

  class PartitioningAlgorithm {
  protected:
    const IPUAlgoConfig& config;

    virtual BatchMapping& createBatch() = 0;
  public:
    std::vector<BatchMapping> mappings;

    PartitioningAlgorithm(const IPUAlgoConfig& config) : config(config), mappings{BatchMapping(config)} {
    }
    virtual bool addComparison(const Comparison& cmpData, BatchMapping& curMapping) = 0;

    std::vector<BatchMapping> getMappings() {
      return std::move(mappings);
    }

    bool addComparison(const Comparison& cmpData) {
      BatchMapping* curMapping = &mappings.back();
      if (!addComparison(cmpData, *curMapping)) {
        curMapping = &createBatch();
        if (!addComparison(cmpData, *curMapping)) {
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

    bool addComparison(const Comparison& cmpData, BatchMapping& curMapping) override {
      for (; bucketIndex < curMapping.buckets.size(); ++bucketIndex) {
        Bucket& bucket = curMapping.buckets[bucketIndex];
        if (bucket.addComparison(cmpData)) {
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
    virtual bool addComparison(const Comparison& cmpData, BatchMapping& curMapping) override {
      BatchMapping& m = mappings.back();
      std::pop_heap(m.buckets.begin(), m.buckets.end(), cmp);
      Bucket& b = m.buckets.back();
      auto success = b.addComparison(cmpData);
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

  std::vector<BatchMapping> mapBatches(IPUAlgoConfig config, const RawSequences& Seqs, Comparisons& Cmps) {

    PLOGD << "Mapping " << Cmps.size() << " comparisons to batches";

    PLOGD << "Set complexity using " << complexityToConfigString(config.complexityAlgo);
    for (auto&& cmp : Cmps) {
      cmp.complexity = calculateComplexity(cmp, config.complexityAlgo);
    }

    if (config.partitioningSortComparisons) {
      PLOGD << "Sorting comparisons";
      std::sort(Cmps.rbegin(), Cmps.rend());
    }

    PLOGD << "Partitioning " << Cmps.size() << " comparison datas";
    auto Partitioner = createPartitioner(config);
    for (int i = 0; i < Cmps.size(); ++i) {
      Partitioner->addComparison(Cmps[i]);
    }
    PLOGD << "Created " << Partitioner->mappings.size() << " batches";
    return Partitioner->getMappings();
  }

SWMeta ComparisonMapping::createMeta() const {
	return {
		.sizeA = comparison->sizeA,
		.offsetA = (int32_t) offsetA,
    .sizeB = comparison->sizeB,
    .offsetB = (int32_t) offsetB,
	};
}
}}