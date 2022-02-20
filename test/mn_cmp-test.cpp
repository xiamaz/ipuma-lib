#include "gtest/gtest.h"

#include "swatlib/vector.hpp"
#include "ipu_batch_affine.h"

class MN_Partitioning : public ::testing::Test {
protected:
  ipu::IPUAlgoConfig algoconfig = {2, 30, 10, 1000};
  std::vector<std::string> seqs = {"AAAA", "TTTT", "GGGG", "CCCC"};
  ipu::Comparisons comparisons = {{0, 1}};
};

using MN_PartitioningDeath = MN_Partitioning;

TEST_F(MN_Partitioning, checkMappingSimple) {
  // ipu::batchaffine::IPUAlgoConfig algoconfig = {2, 30, 10, 1000};
  // std::vector<std::string> seqs = {"AAAA", "TTTT", "GGGG", "CCCC"};
  ipu::partition::BucketMap map(algoconfig.tilesUsed, algoconfig.maxBatches, algoconfig.bufsize);
  ipu::batchaffine::SWAlgorithm::fillMNBuckets(algoconfig.fillAlgo, map, seqs, comparisons);
  EXPECT_EQ(map.numBuckets, 2);
  EXPECT_EQ(map.buckets[0].cmps.size(), 1);
  EXPECT_EQ(map.buckets[1].cmps.size(), 0);
  EXPECT_EQ(map.buckets[0].cmps[0].comparisonIndex, 0);
  EXPECT_EQ(map.buckets[0].cmps[0].offsetA, 0);
  EXPECT_EQ(map.buckets[0].cmps[0].offsetB, 4);
}

// Compare all sequences in A with all sequences in B
TEST(MN_Test, CompareAll) {
	std::vector<std::string> seqs = {
		"AAAAAA",
		"TTTTTT",
		"TTTTAT",
		"AAAAAA",
	};
  ipu::Comparisons comparisons = {
    {0, 2}, // expected 1
    {0, 3}, // expected 6
    {1, 2}, // expected 4
    {1, 3}, // expected 0
  };
  int numWorkers = 1;
  int numCmps = 30;
  int strlen = 20;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0, .gapExtend = -1, .matchValue = 1, .mismatchValue = -1, .ambiguityValue = -1,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::assembly});
  driver.compare_mn_local(seqs, comparisons);
  auto result = driver.get_result();
  EXPECT_EQ(result.scores[0], 1);
  EXPECT_EQ(result.scores[1], 6);
  EXPECT_EQ(result.scores[2], 4);
  EXPECT_EQ(result.scores[3], 0);
}

TEST(MN_Test, CompareAllRR) {
	std::vector<std::string> seqs = {
		"AAAAAA",
		"TTTTTT",
		"TTTTAT",
		"AAAAAA",
	};
  ipu::Comparisons comparisons = {
    {0, 2}, // expected 1
    {0, 3}, // expected 6
    {1, 2}, // expected 4
    {1, 3}, // expected 0
  };
  int numWorkers = 1;
  int numCmps = 30;
  int strlen = 20;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0, .gapExtend = -1, .matchValue = 1, .mismatchValue = -1, .ambiguityValue = -1,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::assembly, ipu::Algorithm::roundRobin});
  driver.compare_mn_local(seqs, comparisons);
  auto result = driver.get_result();
  EXPECT_EQ(result.scores[0], 1);
  EXPECT_EQ(result.scores[1], 6);
  EXPECT_EQ(result.scores[2], 4);
  EXPECT_EQ(result.scores[3], 0);
}

TEST(MN_Test, CompareAllGreedy) {
	std::vector<std::string> seqs = {
		"AAAAAA",
		"TTTTTT",
		"TTTTAT",
		"AAAAAA",
	};
  ipu::Comparisons comparisons = {
    {0, 2}, // expected 1
    {0, 3}, // expected 6
    {1, 2}, // expected 4
    {1, 3}, // expected 0
  };
  int numWorkers = 1;
  int numCmps = 30;
  int strlen = 20;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0, .gapExtend = -1, .matchValue = 1, .mismatchValue = -1, .ambiguityValue = -1,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::assembly, ipu::Algorithm::greedy});
  driver.compare_mn_local(seqs, comparisons);
  auto result = driver.get_result();
  EXPECT_EQ(result.scores[0], 1);
  EXPECT_EQ(result.scores[1], 6);
  EXPECT_EQ(result.scores[2], 4);
  EXPECT_EQ(result.scores[3], 0);
}