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
    ss << "CM[" << comparison.toString() << ": a(o" << offsetA << ") b(o" << offsetB << ")]";
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

  bool Bucket::addComparison(const MultiComparison& md) {
    size_t newTotalLength = totalSequenceLength + md.totalSeqSize;
    size_t newTotalCmps = cmps.size() * NSEEDS + md.comparisonCount * NSEEDS;
    if (newTotalLength <= sequenceCapacity && newTotalCmps <= comparisonCapacity) {
      std::unordered_map<int32_t, size_t> offsets;
      // add sequences
      for (auto [index, size] : md.seqs) {
        seqs.push_back({
          .index = index,
          .offset = totalSequenceLength,
        });
        offsets[index] = totalSequenceLength;
        totalSequenceLength += size;
        longestLength = std::max(longestLength, (size_t) size);
      }

      // add comparisons
      for (const auto& cmp : md.comparisons) {
        cmps.push_back({
          .offsetA = offsets[cmp.indexA],
          .offsetB = offsets[cmp.indexB],
          .comparison = cmp,
        });
        totalCells += cmp.sizeA * cmp.sizeB * NSEEDS;
      }
      return true;
    }
    return false;
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
        .comparison = d
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

  BatchMapping::BatchMapping(int bucketCount, size_t sequenceCapacity, size_t comparisonCapacity) : maxRemainingSpace(sequenceCapacity) {
    for (int i = 0; i < bucketCount; ++i) {
      buckets.push_back(Bucket(i, sequenceCapacity, comparisonCapacity));
    }
  }

  BatchMapping::BatchMapping(IPUAlgoConfig config) : maxRemainingSpace(config.vertexBufferSize) {
    for (int i = 0; i < config.numVertices; ++i) {
      buckets.push_back(Bucket(i, config.vertexBufferSize, config.maxComparisonsPerVertex));
    }
  }

  bool BatchMapping::operator<(const BatchMapping& other) {
    return maxRemainingSpace < other.maxRemainingSpace;
  }

  template<typename C>
  class PartitioningAlgorithm {
  protected:
    const IPUAlgoConfig& config;

    virtual BatchMapping& createBatch() = 0;

    BatchMapping& getMappingHeap() {
      std::pop_heap(mappings.begin(), mappings.end());
      return mappings.back();
    }

    BatchMapping& createMappingHeap() {
      mappings.push_back(BatchMapping(this->config));
      std::push_heap(mappings.begin(), mappings.end());
      return mappings.back();
    }

    void addMappingHeap() {
      std::push_heap(mappings.begin(), mappings.end());
    }
  public:
    std::vector<BatchMapping> mappings;

    PartitioningAlgorithm(const IPUAlgoConfig& config) : config(config), mappings{BatchMapping(config)} {
    }

    virtual bool addComparison(const C& cmpData, BatchMapping& curMapping) = 0;

    std::vector<BatchMapping> getMappings() {
      return std::move(mappings);
    }

    bool addComparisonHeap(const C& cmpData) {
      BatchMapping* curMapping = &getMappingHeap();
      if (!addComparison(cmpData, *curMapping)) {
        addMappingHeap();
        curMapping = &createBatch();
        if (!addComparison(cmpData, *curMapping)) {
          throw std::runtime_error("Could not insert data into a new batch");
        }
      }
      addMappingHeap();
      return true;
    }

    bool addComparison(const C& cmpData) {
      BatchMapping* curMapping = &mappings.back();
      if (!addComparison(cmpData, *curMapping)) {
        curMapping = &createBatch();
        if (!addComparison(cmpData, *curMapping)) {
          throw std::runtime_error("Could not insert data into a new batch");
        }
      }
      return true;
    }

    bool addComparison(const C& cmpDataFront, const C& cmpDataBack) {
      BatchMapping* curMapping = &mappings.back();
      if (!addComparison(cmpDataFront, *curMapping)) {
        if (!addComparison(cmpDataBack, *curMapping)) {
          curMapping = &createBatch();
          if (!addComparison(cmpDataFront, *curMapping)) {
            throw std::runtime_error("Could not insert data into a new batch");
          }
          return true;
        }
        return false;
      }
      return true;
    }
  };

  template<typename C>
  class FillFirstPartitioning : public PartitioningAlgorithm<C> {
    using PartitioningAlgorithm<C>::PartitioningAlgorithm;

    int bucketIndex = 0;

    BatchMapping& createBatch() override {
      this->mappings.push_back(BatchMapping(this->config));
      bucketIndex = 0;
      return this->mappings.back();
    }

    bool addComparison(const C& cmpData, BatchMapping& curMapping) override {
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

  template<typename C>
  class GreedyPartitioning : public PartitioningAlgorithm<C> {
    virtual bool addComparison(const C& cmpData, BatchMapping& curMapping) override {
      // BatchMapping& m = this->mappings.back();
      std::pop_heap(curMapping.buckets.begin(), curMapping.buckets.end(), cmp);
      Bucket& b = curMapping.buckets.back();
      auto success = b.addComparison(cmpData);
      std::push_heap(curMapping.buckets.begin(), curMapping.buckets.end(), cmp);
      return success;
    }

    BatchMapping& createBatch() override {
      this->mappings.push_back(BatchMapping(this->config));
      BatchMapping& m = this->mappings.back();
      std::make_heap(m.buckets.begin(), m.buckets.end(), cmp);
      return m;
    }

  public:
    BucketCompare cmp;
    GreedyPartitioning(const IPUAlgoConfig& config) : PartitioningAlgorithm<C>(config), cmp{} {
      BatchMapping& m = this->mappings.back();
      std::make_heap(m.buckets.begin(), m.buckets.end(), BucketCompare());
    }
  };

  template<typename C>
  std::unique_ptr<PartitioningAlgorithm<C>> createPartitioner(const IPUAlgoConfig& config) {
    switch (config.fillAlgo) {
    case Algorithm::fillFirst:
      return std::unique_ptr<PartitioningAlgorithm<C>>(new FillFirstPartitioning<C>(config));
      break;
    case Algorithm::roundRobin:
      throw std::runtime_error("Not implemented");
      break;
    case Algorithm::greedy:
      return std::unique_ptr<PartitioningAlgorithm<C>>(new GreedyPartitioning<C>(config));
      break;
    }
    throw std::runtime_error("Not any of supported algorithm");
  }

  template<typename C>
  std::vector<BatchMapping> mapBatches(IPUAlgoConfig config, const RawSequences& Seqs, std::vector<C>& Cmps) {
    addComplexity(Cmps, config.complexityAlgo);

    if (config.partitioningSortComparisons) {
      PLOGD << "Sorting comparisons";
      std::sort(Cmps.rbegin(), Cmps.rend());
    }

    PLOGD << "Partitioning " << Cmps.size() << " comparison datas";
    auto Partitioner = createPartitioner<C>(config);
    switch (config.partitionadd) {
    case PartitionAdd::alternating:
      for (int i = 0, j = Cmps.size() - 1; i <= j;) {
        bool consumedFront = Partitioner->addComparison(Cmps[i], Cmps[j]);
        if (consumedFront) i++; else j--;
      }
      break;
    case PartitionAdd::sequential:
      for (int i = 0; i < Cmps.size(); ++i) {
        bool consumedFront = Partitioner->addComparison(Cmps[i]);
      }
      break;
    case PartitionAdd::heap:
      for (int i = 0; i < Cmps.size(); ++i) {
        bool consumedFront = Partitioner->addComparisonHeap(Cmps[i]);
      }
      break;
    }
    PLOGD << "Created " << Partitioner->mappings.size() << " batches";
    return Partitioner->getMappings();
  }

  template
  std::vector<BatchMapping> mapBatches<Comparison>(IPUAlgoConfig config, const RawSequences& Seqs, Comparisons& Cmps);

  template
  std::vector<BatchMapping> mapBatches<MultiComparison>(IPUAlgoConfig config, const RawSequences& Seqs, MultiComparisons& Cmps);

SWMeta ComparisonMapping::createMeta() const {
	return {
		.sizeA = comparison.sizeA,
		.offsetA = (int32_t) offsetA,
    .sizeB = comparison.sizeB,
    .offsetB = (int32_t) offsetB,
	};
}
}}