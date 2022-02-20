#include "gtest/gtest.h"

#include "ipu_batch_affine.h"

class PartitioningTest : public ::testing::Test {
protected:
  void SetUp() override {
    map = ipu::partition::BucketMap(tilesUsed, maxBatches, bufsize);
  }

  void genFull() {
    for (int i = 0; i < tilesUsed * maxBatches; ++i) {
      a.push_back(std::string("A", maxAB));
      b.push_back(std::string("A", maxAB));
    }
  }

  void checkBucketIndices() {
    EXPECT_EQ(map.buckets.size(), tilesUsed);
    for (int i = 0; i < tilesUsed; ++i) {
      const auto& bucket = map.buckets[i];
      EXPECT_EQ(bucket.bucketIndex, i);
    }
  }

  int tilesUsed = 10;
  int maxBatches = 2;
  int maxAB = 10;
  int bufsize = maxBatches * maxAB * 2;

  std::vector<std::string> a, b;

  ipu::partition::BucketMap map;
};

TEST_F(PartitioningTest, fullFillFirst) {
  genFull();
  ipu::batchaffine::SWAlgorithm::fillBuckets(ipu::Algorithm::fillFirst, map, a, b, 0);

  checkBucketIndices();
  for (int i = 0; i < tilesUsed; ++i) {
    const auto& bucket = map.buckets[i];
    EXPECT_EQ(bucket.maxLen, maxAB);
    EXPECT_EQ(bucket.seqSize, maxAB * maxBatches * 2);
  }
}

TEST_F(PartitioningTest, fullRoundRobin) {
  genFull();
  ipu::batchaffine::SWAlgorithm::fillBuckets(ipu::Algorithm::roundRobin, map, a, b, 0);

  checkBucketIndices();
  for (int i = 0; i < tilesUsed; ++i) {
    const auto& bucket = map.buckets[i];
    EXPECT_EQ(bucket.maxLen, maxAB);
    EXPECT_EQ(bucket.seqSize, maxAB * maxBatches * 2);
  }
}

TEST_F(PartitioningTest, fullGreedy) {
  genFull();
  ipu::batchaffine::SWAlgorithm::fillBuckets(ipu::Algorithm::greedy, map, a, b, 0);

  checkBucketIndices();
  for (int i = 0; i < tilesUsed; ++i) {
    const auto& bucket = map.buckets[i];
    EXPECT_EQ(bucket.maxLen, maxAB);
    EXPECT_EQ(bucket.seqSize, maxAB * maxBatches * 2);
  }
}

TEST_F(PartitioningTest, partialFillFirst) {
  a = {"AAAAAA", "AA"};
  b = {"AAAAAA", "AAA"};
  ipu::batchaffine::SWAlgorithm::fillBuckets(ipu::Algorithm::fillFirst, map, a, b, 0);

  checkBucketIndices();
  for (int i = 0; i < tilesUsed; ++i) {
    const auto& bucket = map.buckets[i];
    if (i > 0) {
      EXPECT_EQ(bucket.maxLen, 0);
      EXPECT_EQ(bucket.seqSize, 0);
    } else {
      EXPECT_EQ(bucket.maxLen, 6);
      EXPECT_EQ(bucket.seqSize, 17);
    }
  }
}

TEST_F(PartitioningTest, partialRoundRobin) {
  a = {"AAAAAA", "AA"};
  b = {"AAAAAA", "AAA"};
  ipu::batchaffine::SWAlgorithm::fillBuckets(ipu::Algorithm::roundRobin, map, a, b, 0);

  checkBucketIndices();
  for (int i = 0; i < tilesUsed; ++i) {
    const auto& bucket = map.buckets[i];
    if (i > 1) {
      EXPECT_EQ(bucket.maxLen, 0);
      EXPECT_EQ(bucket.seqSize, 0);
    } else if (i == 0) {
      EXPECT_EQ(bucket.maxLen, 6);
      EXPECT_EQ(bucket.seqSize, 12);
    } else if (i == 1) {
      EXPECT_EQ(bucket.maxLen, 3);
      EXPECT_EQ(bucket.seqSize, 5);
    }
  }
}

TEST(PartitionMNTest, compression) {
  ipu::RawSequences Seqs = {"AAAA", "AAAA", "AAAA"};
  ipu::Comparisons Cmps = {{0, 1}, {0, 2}, {0, 2}};
  ipu::partition::BucketMap map(1, 5, 100);
  ipu::batchaffine::SWAlgorithm::fillMNBuckets(ipu::Algorithm::fillFirst, map, Seqs, Cmps);
  EXPECT_EQ(map.buckets[0].seqSize, 12);
}

TEST(PartitionMNTest, compressionRR) {
  ipu::RawSequences Seqs = {"AAAA", "AAAA", "AAAA"};
  ipu::Comparisons Cmps = {{0, 1}, {0, 2}, {0, 2}};
  ipu::partition::BucketMap map(1, 5, 100);
  ipu::batchaffine::SWAlgorithm::fillMNBuckets(ipu::Algorithm::roundRobin, map, Seqs, Cmps);
  EXPECT_EQ(map.buckets[0].seqSize, 12);
}

TEST(PartitionMNTest, compressionGreedy) {
  ipu::RawSequences Seqs = {"AAAA", "AAAA", "AAAA"};
  ipu::Comparisons Cmps = {{0, 1}, {0, 2}, {0, 2}};
  ipu::partition::BucketMap map(1, 5, 100);
  ipu::batchaffine::SWAlgorithm::fillMNBuckets(ipu::Algorithm::greedy, map, Seqs, Cmps);
  EXPECT_EQ(map.buckets[0].seqSize, 12);
}