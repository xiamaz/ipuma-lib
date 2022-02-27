#include "gtest/gtest.h"

#include "ipu_batch_affine.h"
#include "swatlib/vector.hpp"

using std::max;
using std::min;
using std::string;
using std::vector;

TEST(IPUDev, DISABLED_TestMultiStriping) {
  int numWorkers = 1;
  int numCmps = 6;
  int strlen = 6;
  int bufsize = 6;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0,
    .gapExtend = -ALN_GAP_EXTENDING_COST,
    .matchValue = ALN_MATCH_SCORE,
    .mismatchValue = -ALN_MISMATCH_COST,
    .ambiguityValue = -ALN_AMBIGUITY_COST,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::multiasm});

  vector<string> queries, refs;
  for (int i = 0; i < 1; ++i) {
    queries.push_back(std::string(strlen, 'A'));
    refs.push_back(std::string(strlen, 'A'));
    // queries.push_back("AAAAAAAAAA");
    // refs.push_back("AAAAAAAAAA");
  }
  // std::cout << queries[0] << "\n";
  // refs[1] = "TTAAAA";
  // refs[4] = "TTTTTT";

  driver.compare_local(queries, refs, false);
  auto aln_results = driver.get_result();
  // std::cout << driver.printTensors();
  std::cout << "Aln results:\n";
  std::cout << swatlib::printVectorD(aln_results.scores) << "\n";
  std::cout << swatlib::printVectorD(aln_results.a_range_result) << "\n";
}

TEST(IPUDev, DISABLED_TestStripingAsm) {
  int numWorkers = 1;
  int numCmps = 6;
  int strlen = 6;
  int bufsize = 6;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0,
    .gapExtend = -ALN_GAP_EXTENDING_COST,
    .matchValue = ALN_MATCH_SCORE,
    .mismatchValue = -ALN_MISMATCH_COST,
    .ambiguityValue = -ALN_AMBIGUITY_COST,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::multiasm});

  vector<string> queries, refs;
  for (int i = 0; i < 1; ++i) {
    queries.push_back(std::string(strlen, 'A'));
    refs.push_back(std::string(strlen, 'A'));
    // queries.push_back("AAAAAAAAAA");
    // refs.push_back("AAAAAAAAAA");
  }
  // std::cout << queries[0] << "\n";
  // refs[1] = "TTAAAA";
  // refs[4] = "TTTTTT";

  driver.compare_local(queries, refs, false);
  auto aln_results = driver.get_result();
  // std::cout << driver.printTensors();
  std::cout << "Aln results:\n";
  std::cout << swatlib::printVectorD(aln_results.scores) << "\n";
  std::cout << swatlib::printVectorD(aln_results.a_range_result) << "\n";
}

TEST(IPUDev, DISABLED_MultiVertexSeparate) {
  int numWorkers = 1;
  int numCmps = 30;
  int strlen = 20;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0,
    .gapExtend = -ALN_GAP_EXTENDING_COST,
    .matchValue = ALN_MATCH_SCORE,
    .mismatchValue = -ALN_MISMATCH_COST,
    .ambiguityValue = -ALN_AMBIGUITY_COST,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::multiasm});

  vector<string> queries, refs;
  for (int i = 0; i < 6; ++i) {
    queries.push_back("AAAAAA");
    refs.push_back("AAAAAA");
  }
  refs[1] = "TTAAAA";
  refs[4] = "TTTTTT";
  driver.compare_local(queries, refs);
  auto aln_results = driver.get_result();
  std::cout << "Aln results:\n";
  std::cout << swatlib::printVector(aln_results.scores) << "\n";
}