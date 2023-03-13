#include "gtest/gtest.h"
#include <plog/Log.h>

#include "ipu_batch_affine.h"

class PartitioningTest : public ::testing::Test {
protected:
  ipu::RawSequences Seqs;
  ipu::Comparisons Cmps;

  ipu::IPUAlgoConfig config;
};

/**
 * @brief Test perfectly filling a batch
 * 
 */
TEST_F(PartitioningTest, fullFillFirst) {
  int sequence_length = 100;
  config = {
    3 /*Tiles used*/,
    1000 /*maxSequenceLength*/,
    200000,
    sequence_length * 2 * 10,
    ipu::VertexType::multixdrop,
    ipu::Algorithm::fillFirst,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  for (int i = 0; i < numberCmps; ++i) {
    Seqs.push_back(std::string(sequence_length, 'A'));
    Seqs.push_back(std::string(sequence_length, 'B'));
    Cmps.push_back({.indexA = i * 2, .indexB = i * 2 + 1});
  }

  auto batches = ipu::partition::mapBatches(config, Seqs, Cmps);

  ASSERT_EQ(batches.size(), 1);
}

TEST_F(PartitioningTest, fullGreedy) {
  int sequence_length = 100;
  config = {
    3 /*Tiles used*/,
    1000 /*maxSequenceLength*/,
    200000,
    sequence_length * 2 * 10,
    ipu::VertexType::multixdrop,
    ipu::Algorithm::greedy,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  for (int i = 0; i < numberCmps; ++i) {
    Seqs.push_back(std::string(sequence_length, 'A'));
    Seqs.push_back(std::string(sequence_length, 'B'));
    Cmps.push_back({.indexA = i * 2, .indexB = i * 2 + 1});
  }

  auto batches = ipu::partition::mapBatches(config, Seqs, Cmps);

  ASSERT_EQ(batches.size(), 1);
}


/**
 * @brief Test perfectly filling two batches
 * 
 */
TEST_F(PartitioningTest, MultiBatchSizedCount) {
  int sequence_length = 100;
  int batches_count = 4;
  config = {
    3 /*Tiles used*/,
    1000 /*maxSequenceLength*/,
    200000,
    sequence_length * 2 * 10,
    ipu::VertexType::multixdrop,
    ipu::Algorithm::fillFirst,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  for (int i = 0; i < numberCmps * batches_count; ++i) {
    Seqs.push_back(std::string(sequence_length, 'A'));
    Seqs.push_back(std::string(sequence_length, 'B'));
    Cmps.push_back({.indexA = i * 2, .sizeA = (int32_t) Seqs[i*2].size(), .indexB = i * 2 + 1, .sizeB = (int32_t) Seqs[i*2+1].size()});
  }


  auto batches = ipu::partition::mapBatches(config, Seqs, Cmps);
  for (auto && batch : batches) {
  size_t bufsum = 0;
    for (auto && bucket : batch.buckets) {
      for (auto && cmp : bucket.cmps) {
        bufsum += cmp.comparison->sizeA;
        bufsum += cmp.comparison->sizeB;
      }
    }
    ASSERT_EQ(bufsum,  2 * numberCmps * sequence_length);
  }
  ASSERT_EQ(batches.size(), batches_count);
}