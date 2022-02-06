#include "gtest/gtest.h"

#include "ipu_batch_affine.h"


// Compare all sequences in A with all sequences in B
TEST(MN_Test, CompareAll) {
	std::vector<std::string> A = {
		"AAAAAA",
		"AAAAAA",
	};
	std::vector<std::string> B = {
		"AAAAAA",
		"TTTTTT",
	};
  std::vector<int> comparisons = {
    0, 0,
    0, 1,
    1, 0,
    1, 1,
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
  driver.compare_mn_local(A, B, comparisons);
  auto result = driver.get_result();
}