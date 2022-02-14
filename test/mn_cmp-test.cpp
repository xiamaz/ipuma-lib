#include "gtest/gtest.h"

#include "swatlib/vector.hpp"
#include "ipu_batch_affine.h"

class MN_Partitioning : public ::testing::Test {
protected:
  ipu::batchaffine::IPUAlgoConfig algoconfig = {2, 30, 10, 1000};
  std::vector<std::string> seqs = {"AAAA", "TTTT", "GGGG", "CCCC"};
  std::vector<int> comparisons = {0, 1};
};

using MN_PartitioningDeath = MN_Partitioning;

TEST_F(MN_Partitioning, checkMappingSimple) {
  // ipu::batchaffine::IPUAlgoConfig algoconfig = {2, 30, 10, 1000};
  // std::vector<std::string> seqs = {"AAAA", "TTTT", "GGGG", "CCCC"};
  auto buckets = ipu::batchaffine::SWAlgorithm::fillMNBuckets(algoconfig, seqs, comparisons);
  EXPECT_EQ(buckets.size(), 2);
  EXPECT_EQ(buckets[0].comparisons.size(), 1);
  EXPECT_EQ(buckets[1].comparisons.size(), 0);
  EXPECT_EQ(buckets[0].comparisons[0].cmpIndex, 0);
  EXPECT_EQ(buckets[0].comparisons[0].indexA, 0);
  EXPECT_EQ(buckets[0].comparisons[0].indexB, 1);
}

TEST_F(MN_Partitioning, checkInvalidComparisonLength) {
  comparisons = {0, 1, 2}; // comparisons needs to be always a multiple of 2
  EXPECT_THROW({
    ipu::batchaffine::SWAlgorithm::fillMNBuckets(algoconfig, seqs, comparisons);
  }, std::logic_error);
}

// Compare all sequences in A with all sequences in B
TEST(MN_Test, CompareAll) {
	std::vector<std::string> seqs = {
		"AAAAAA",
		"TTTTTT",
		"TTTTAT",
		"AAAAAA",
	};
  std::vector<int> comparisons = {
    0, 2, // expected 1
    0, 3, // expected 6
    1, 2, // expected 4
    1, 3, // expected 0
  };
  int numWorkers = 1;
  int numCmps = 30;
  int strlen = 20;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0, .gapExtend = -1, .matchValue = 1, .mismatchValue = -1, .ambiguityValue = -1,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::batchaffine::VertexType::assembly});
  driver.compare_mn_local(seqs, comparisons);
  auto result = driver.get_result();
  std::cout << swatlib::printVector(result.scores) << "\n";
}