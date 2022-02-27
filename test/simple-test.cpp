#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "ipu_batch_affine.h"
#include "swatlib/vector.hpp"

using std::max;
using std::min;
using std::string;
using std::vector;

struct SequenceTestData {
  string ref;
  string query;
  int qstart;
  int qend;
  int rstart;
  int rend;
  int score;
};

const vector<SequenceTestData> simpleData = {
  {"acgt", "acgt", 0, 3, 0, 3, 4},
  {"ACGT", "ACGT", 0, 3, 0, 3, 4},
  {"AACGT", "ACGT", 0, 3, 1, 4, 4},
  {"ACGTT", "ACGT", 0, 3, 0, 3, 4},
  {"ACGT", "TACGT", 1, 4, 0, 3, 4},
  {"ACGT", "TTACGT", 2, 5, 0, 3, 4},
  {"ACGT", "ACGTTT", 0, 3, 0, 3, 4},
  {"ACGT", "TACGTT", 1, 4, 0, 3, 4},
  {"ACGT", "TTACGTT", 2, 5, 0, 3, 4},
  {"ACGT", "TACGTTT", 1, 4, 0, 3, 4},
  {"ACGT", "TTACGTTT", 2, 5, 0, 3, 4},
  {"ACGT", "TTACGTTT", 2, 5, 0, 3, 4},
  {"AAAATTTTCCCCGGGG", "AAAATTTTCCCCGGGG", 0, 15, 0, 15, 16},
  {"AAAATTTTCCCCGGGG", "AAAATTTTACCCGGGG", 0, 15, 0, 15, 14},
  {"AAAATTTTCCCCGGGG", "AAAATTTTACCCCGGGG", 0, 16, 0, 15, 15},
  {"AAAATTTTCCCCGGGG", "AAAATTTCCCCGGGG", 0, 14, 0, 15, 14},
  {"AAAATTTTCCCCGGGG", "GCTAGCTAGCTAGCTA", 3, 3, 0, 0, 1},
  {"AAAATTTTCCCCGGGG", "GCTAAAATTTTCCCCGGGG", 3, 18, 0, 15, 16},
  {"AAAATTTTCCCCGGGG", "AAAATTTTCCCCGGGGACT", 0, 15, 0, 15, 16},
};

class SimpleCorrectnessTest : public ::testing::Test {
protected:
  void SetUp() override {
    std::transform(testData.begin(), 
                   testData.end(), 
                   std::back_inserter(queries), 
                  [](const SequenceTestData& t) { return t.query; });
    std::transform(testData.begin(), 
                   testData.end(), 
                   std::back_inserter(refs), 
                  [](const SequenceTestData& t) { return t.ref; });
  }

  void checkResults(ipu::batchaffine::BlockAlignmentResults results) {
    for (int i = 0; i < testData.size(); ++i) {
      const auto& [ref, query, qstart, qend, rstart, rend, score] = testData[i];

      auto uda = results.a_range_result[i];
      int16_t query_begin = uda & 0xffff;
      int16_t query_end = uda >> 16;

      auto udb = results.b_range_result[i];
      int16_t ref_begin = udb & 0xffff;
      int16_t ref_end = udb >> 16;

      EXPECT_EQ(results.scores[i], score) << i << ": Alignment score does not match expected value";
      EXPECT_EQ(ref_begin, rstart) << i << ": Alignment reference start does not match expected value";
      EXPECT_EQ(ref_end, rend) << i << ": Alignment reference end does not match expected value";
      EXPECT_EQ(query_begin, qstart) << i << ": Alignment query start does not match expected value";
      EXPECT_EQ(query_end, qend) << i << ": Alignment query end does not match expected value";
    }
  }

  const vector<SequenceTestData>& testData = simpleData;
  vector<string> queries;
  vector<string> refs;
};

TEST_F(SimpleCorrectnessTest, UseAssemblyVertex) {
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
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::assembly});

  driver.compare_local(queries, refs);
  auto aln_results = driver.get_result();
  checkResults(aln_results);
}

TEST_F(SimpleCorrectnessTest, UseCppVertex) {
  auto driver = ipu::batchaffine::SWAlgorithm({}, {
    .tilesUsed = 2,
    .maxAB = 300,
    .maxBatches = 20,
    .bufsize = 3000,
    .vtype = ipu::VertexType::cpp,
    .fillAlgo = ipu::Algorithm::roundRobin
  });

  driver.compare_local(queries, refs);
  auto aln_results = driver.get_result();
  checkResults(aln_results);
}

TEST_F(SimpleCorrectnessTest, UseCppMultiVertex) {
  auto driver = ipu::batchaffine::SWAlgorithm({}, {
    .tilesUsed = 2,
    .maxAB = 300,
    .maxBatches = 20,
    .bufsize = 3000,
    .vtype = ipu::VertexType::multi,
    .fillAlgo = ipu::Algorithm::roundRobin
  });

  driver.compare_local(queries, refs);
  auto aln_results = driver.get_result();
  checkResults(aln_results);
}

TEST_F(SimpleCorrectnessTest, UseAsmMultiVertex) {
  auto driver = ipu::batchaffine::SWAlgorithm({}, {
    .tilesUsed = 2,
    .maxAB = 300,
    .maxBatches = 20,
    .bufsize = 3000,
    .vtype = ipu::VertexType::multiasm,
    .fillAlgo = ipu::Algorithm::roundRobin
  });

  driver.compare_local(queries, refs);
  auto aln_results = driver.get_result();
  checkResults(aln_results);
}

TEST_F(SimpleCorrectnessTest, useMNComparisons) {
  auto driver = ipu::batchaffine::SWAlgorithm({}, {
    .tilesUsed = 4,
    .maxAB = 300,
    .maxBatches = 5,
    .bufsize = 3000,
    .vtype = ipu::VertexType::cpp,
    .fillAlgo = ipu::Algorithm::roundRobin
  });

  std::vector<std::string> seqs;
  ipu::Comparisons comparisons;
  for (int i = 0; i < queries.size(); ++i) {
    seqs.push_back(queries[i]);
    seqs.push_back(refs[i]);
    comparisons.push_back({2 * i, 2 * i + 1});
  }

  driver.compare_mn_local(seqs, comparisons);
  auto aln_results = driver.get_result();
  checkResults(aln_results);
}

TEST_F(SimpleCorrectnessTest, prepareTest) {
  int numWorkers = 1;
  int numCmps = 30;
  int strlen = 20;
  int bufsize = 1000;
  ipu::SWConfig swconfig = {
    .gapInit = 0,
    .gapExtend = -ALN_GAP_EXTENDING_COST,
    .matchValue = ALN_MATCH_SCORE,
    .mismatchValue = -ALN_MISMATCH_COST,
    .ambiguityValue = -ALN_AMBIGUITY_COST,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  };
  ipu::IPUAlgoConfig algoconfig = {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::assembly};

  std::vector<int32_t> inputs(algoconfig.getInputBufferSize32b());
  std::vector<int> mapping(queries.size(), 0);

  ipu::batchaffine::SWAlgorithm::prepare_remote(swconfig, algoconfig, queries, refs, &*inputs.begin(), &*inputs.end(), mapping.data());
}

TEST(PrepareTest, simple) {
  int tilesUsed = 10;
  int maxBatches = 2;
  int maxAB = 10;
  int bufsize = maxBatches * maxAB * 2;
  ipu::SWConfig swconfig = {};
  ipu::IPUAlgoConfig config = {
    .tilesUsed = tilesUsed,
    .maxAB = maxAB,
    .maxBatches = maxBatches,
    .bufsize = bufsize,
    .vtype = ipu::VertexType::cpp,
    .fillAlgo = ipu::Algorithm::fillFirst
  };

  std::vector<std::string> A, B;
  for (int i = 0; i < tilesUsed * maxBatches; ++i) {
    A.push_back(std::string(maxAB, 'A'));
    B.push_back(std::string(maxAB, 'A'));
  }

  std::vector<int32_t> inputs(config.getInputBufferSize32b() + 2);
  std::vector<int> mapping(A.size(), 0);

  inputs[0] = 0xDEADBEEF;
  inputs[inputs.size() - 1] = 0xDEADBEEF;

  ipu::batchaffine::SWAlgorithm::prepare_remote(swconfig, config, A, B, &*inputs.begin() + 1, &*inputs.end() - 1, mapping.data());
  std::vector<int32_t> slice_begin(inputs.begin(), inputs.begin() + 10);
  std::vector<int32_t> slice_end(inputs.end() - 10, inputs.end());
  EXPECT_EQ(*inputs.begin(), 0xDEADBEEF) << "Start overwritten: " << swatlib::printVector(slice_begin);
  EXPECT_EQ(*(inputs.end() - 1), 0xDEADBEEF) << "End overwritten: " << swatlib::printVector(slice_end);
}

TEST(LargerTest, simple) {
  int tilesUsed = 10;
  int maxBatches = 2;
  int maxAB = 10;
  int bufsize = maxBatches * maxAB * 2;
  ipu::SWConfig swconfig = {};
  ipu::IPUAlgoConfig config = {
    .tilesUsed = tilesUsed,
    .maxAB = maxAB,
    .maxBatches = maxBatches,
    .bufsize = bufsize,
    .vtype = ipu::VertexType::assembly,
    .fillAlgo = ipu::Algorithm::fillFirst
  };

  std::vector<std::string> A, B;
  for (int i = 0; i < 2; ++i) {
    A.push_back(std::string(maxAB, 'A'));
    B.push_back("TTTT" + std::string(maxAB, 'A'));
  }
  auto driver = ipu::batchaffine::SWAlgorithm(swconfig, config);
  driver.compare_local(A, B);
  auto result = driver.get_result();
  std::cout << A[0].size() << ": " << swatlib::printVector(result.scores) << "\n";
}