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
    Seqs.push_back(std::string("A", sequence_length));
    Seqs.push_back(std::string("B", sequence_length));
    Cmps.push_back({.indexA = i * 2, .indexB = i * 2 + 1});
  }

  auto batches = ipu::partition::mapBatches(config, Seqs, Cmps);

  ASSERT_EQ(batches.size(), 1);
}

TEST_F(PartitioningTest, fullGreedy) {
  size_t sequence_length = 100;
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
    Seqs.push_back(std::string("A", sequence_length));
    Seqs.push_back(std::string("B", sequence_length));
    Cmps.push_back({.indexA = i * 2, .indexB = i * 2 + 1});
  }

  auto batches = ipu::partition::mapBatches(config, Seqs, Cmps);

  ASSERT_EQ(batches.size(), 1);
}