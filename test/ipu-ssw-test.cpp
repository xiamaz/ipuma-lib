#include <sstream>
#include <string>
#include <vector>

#include "ssw.hpp"
#include "gtest/gtest.h"

#include "ipu_batch_affine.h"
#include "swatlib/vector.hpp"

using std::max;
using std::min;
using std::string;
using std::vector;

using namespace StripedSmithWaterman;

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

void test_aligns_ipu(vector<Alignment> &alns, vector<string> query, vector<string> ref, ipu::batchaffine::SWAlgorithm &algo) {
  alns.reserve(query.size());

  algo.compare_local(query, ref);
  auto aln_results = algo.get_result();

  for (int i = 0; i < query.size(); i++) {
    alns.push_back({});
    Alignment &aln = alns[i];

    auto uda = aln_results.a_range_result[i];
    int16_t query_begin = uda & 0xffff;
    int16_t query_end = uda >> 16;

    auto udb = aln_results.b_range_result[i];
    int16_t ref_begin = udb & 0xffff;
    int16_t ref_end = udb >> 16;

    aln.query_end = query_end;
    aln.query_begin = query_begin;
    aln.ref_begin = ref_begin;
    aln.ref_end = ref_end;
    aln.sw_score = aln_results.scores[i];
    aln.sw_score_next_best = 0;
    aln.mismatches = 0;
  }
}

inline string aln2string(Alignment &aln) {
  std::stringstream ss;
  ss << "score=" << aln.sw_score << " score2=" << aln.sw_score_next_best;
  ss << " rbegin=" << aln.ref_begin << " rend=" << aln.ref_end;
  ss << " qbegin=" << aln.query_begin << " qend=" << aln.query_end;
  ss << " rend2=" << aln.ref_end_next_best << " mismatches=" << aln.mismatches;
  ss << " cigarstr=" << aln.cigar_string;
  return ss.str();
}

class ParityTest : public ::testing::Test {
protected:
  void SetUp() override {
  queries = {
    "AATGAGAATGATGTCGTTCGAAATTTGACCAGTCAAACCGCGGGCAATAAGGTCTTCGTTCAGGGCATAGACCTTAATGGGGGCATTACGCAGACTTTCA",
    "ATCTGGCAGGTAAAGATGAGCTCAACAAAGTGATCCAGCATTTTGGCAAAGGAGGCTTTGATGTGATTACTCGCGGTCAGGTGCCACCTAACCCGTCTGA",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGGTATTGGGAGCACG",
    "CCCCGCACCCGCAAGCCGCCGAGAAAAAAAGGATGAGGGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGAACA",
    "CGTGCTCCCAATACCCGCCAGGAATTAGCGCGCTTACGCTCAGAGCTGGCGGAAAAAGAGAGCCAGCAACTGGATGCGTATTTGCAGGCCTTGCGTAATC",
    "CCCCGCACCCGCAAGCCGCCGAGAAAAAAAGGATGAGGGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGAACA",
    "ATATCATCACTCCGATGGACGTTTCGATTATCGGCTGCGTGGGGAATGGCCCAGGTGAGGCGCTGGTTTCTACACTCGGCGTCACCGGCGGGACGAAACA",
    "CCGGCGGTGGGCGCGTCCGCCAGTGCCGGCGCGAGCAGGACGGCGTAGAGCCGGATGACGTGATCCTGCCGCCGGAGAAGGCAAAAGCGAAAGGGCTGAC",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGGTATTGGGAGCACG",
    "TCCGGCTGGCAGAACTTGACCAGTGCCGATAAAGAAAAGATGTGGCAGGCGCAGGGGCGAATGCTCACCGCACAAAGCCTGAAGATTAATGCGTTGCTGC",
    "CCGCCCCCGCGCACACGGTGCGGCCTGTCCCGCGTATACTCGACCAGCGGCGTCCCGCCCAGCTTCATTCCCGCCAGGTAACCGCTGCCATACGTCAGCC",
    "AAACCATTGCGAGGAAGTGGTTCTACTTGCTGTCGCCGCGGGAGAACAGGTTGTGCCTGCCTCTGAACTTGCTGCCGCCATGAAGCAGATTAAAGAACTC",
    "CCTTCCCCCTAACTTTCCGCCCGCCATGAAGCAGATAAAAGAACTCCAGCGCCTGCTCGGAAAGAAAACGATGGAAAATGAACTCCTCAAAGAAGCCGTT",
    "AGATGTGCCGGTCATTAAGCATAAAGCCGATGGTTTCTCCCCGCACTTGCCGCCAGTGACGCCACGGCCAGTCAGAGAAGATCATAACAACCGCTCCAGT",
    "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCCGTAGTTCTTCACCCGTTTTGAATTCGGGTTTGTAT",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCACTTGATGGCAGATAACGCGCC",
  };
  refs = {
    "AATGAGAATGATGTCNTTCNAAATTTGACCAGTCAAACCGCGGGCAATAAGGTCTTCGTTCAGGGCATAGACCTTAATGGGGGCATTACGCAGACTTTCA",
    "ATCTGGCAGGTAAAGATGAGCTCAACAAAGTGATCCAGCATTTTGGCAAAGGAGGCTTTGATGTGATTACTCGCGGTCAGGTGCCACCTAANNCGTCTGA",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGTTATTGGCAGACAG",
    "GCACCGTCCAGCCAACCGCCGAGAAGAAAAGAATGAGTGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGACAC",
    "CTGTCTGCCAATAACCGCCAGGAATTAGCGCGCTTACGCTCAGAGCTGGCGGAAAAAGAGAGCCAGCAACTGGATGCGTATTTGCAGGCCTTGCGTAATC",
    "GCACCGTCCAGCCAACCGCCGAGAAGAAAAGAATGAGTGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGACAC",
    "TGGTTTCTACACTCGGCGTCACCGGCGGCAACAAGAA",
    "AGCGCCGGGCGCGCTTCCGCCAGTGCCTGCGCGAGCAGGACGGCGTAGAGCCGGATGACGTGATCCTGCCGCCGGAGAAGGCAAAAGCGAAAGGGCTGAC",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGTTATTGGCAGACAG",
    "TTCGCCGCGCAGAACCTGACCAGTGCCGATAACGAAAAGATGTGGCAGGCGCAGGGGCGAATGCTCACCGCACAAAGCCTGAAGATTAATGCGTTGCTGC",
    "CTGCGCACCGTCTCACGGTGCAGCCTGTCCCGCGTATACTCGACCAGCGGCGTCCCGCCCAGCTTCATTCCCGCCAGGTAACCGCTGCCATACGTCAGC",
    "TAAGCAATACCAGGAAGGAAGTCTTACTGCTGTCGCCGCCGGAGAACAGGTTGTTCCTGCCTCTGAACTTGCTGCCGCCATGAAGCAGATTAAAGAACTC",
    "TCCTGCCTCTGAACTTGCTGCCGCCATGAAGCAGATTAAAGAACTCCAGCGCCTGCTCGGCAAGAAAACGATGGAAAATGAACTCCTCAAAGAAGCCGTT",
    "TGAGTTGCTCGTCATTAAGACGTAAGGCGATGGTTTCTCCCCGCACTTGCCGCCAGTG",
    "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCAGTAGTTCTTCACCCGTTTTGAATGCGGGTTTGATA",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCAACTGCTGGCAGATAACCCCGC",
  };
  }

  vector<string> queries, refs;
  
  void checkResults(const vector<Alignment>& alns_ipu) {

    Aligner ssw_aligner(
      ALN_MATCH_SCORE,
      ALN_MISMATCH_COST,
      ALN_GAP_OPENING_COST,
      ALN_GAP_EXTENDING_COST,
      ALN_AMBIGUITY_COST);

    StripedSmithWaterman::Filter ssw_filter;
    vector<Alignment> alns(queries.size());
    for (int i = 0; i < queries.size(); ++i) {
      auto reflen = refs[i].size();
      auto qlen = queries[i].size();

      ssw_aligner.Align(refs[i].data(), reflen, queries[i].data(), qlen, ssw_filter, &alns[i], max((int)(qlen / 2), 15));
      EXPECT_EQ(alns_ipu[i].sw_score, alns[i].sw_score) << i << ": IPU score result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].query_begin, alns[i].query_begin) << i << ": IPU reference start result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].query_end, alns[i].query_end) << i << ": IPU reference end result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].ref_begin, alns[i].ref_begin) << i << ": IPU query start result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].ref_end, alns[i].ref_end) << i << ": IPU query end result does not match CPU SSW";
    }
  }
};

TEST_F(ParityTest, UseAssembly) {
  int numWorkers = 10;
  int numCmps = 10;
  int strlen = 120;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -(ALN_GAP_OPENING_COST-ALN_GAP_EXTENDING_COST),
    .gapExtend = -ALN_GAP_EXTENDING_COST,
    .matchValue = ALN_MATCH_SCORE,
    .mismatchValue = -ALN_MISMATCH_COST,
    .ambiguityValue = -ALN_AMBIGUITY_COST,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::batchaffine::VertexType::assembly, ipu::batchaffine::partition::Algorithm::greedy});

  vector<Alignment> alns_ipu(queries.size());
  test_aligns_ipu(alns_ipu, refs, queries, driver);

  checkResults(alns_ipu);
}

TEST_F(ParityTest, UseMultiAssembly) {
  int numWorkers = 10;
  int numCmps = 10;
  int strlen = 120;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -(ALN_GAP_OPENING_COST-ALN_GAP_EXTENDING_COST),
    .gapExtend = -ALN_GAP_EXTENDING_COST,
    .matchValue = ALN_MATCH_SCORE,
    .mismatchValue = -ALN_MISMATCH_COST,
    .ambiguityValue = -ALN_AMBIGUITY_COST,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::batchaffine::VertexType::multiasm, ipu::batchaffine::partition::Algorithm::greedy});

  vector<Alignment> alns_ipu(queries.size());
  test_aligns_ipu(alns_ipu, refs, queries, driver);

  checkResults(alns_ipu);
}

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
  }, {numWorkers, strlen, numCmps, bufsize, ipu::batchaffine::VertexType::multistripedasm});

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
  }, {numWorkers, strlen, numCmps, bufsize, ipu::batchaffine::VertexType::stripedasm});

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
  }, {numWorkers, strlen, numCmps, bufsize, ipu::batchaffine::VertexType::multiasm});

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
  }, {numWorkers, strlen, numCmps, bufsize, ipu::batchaffine::VertexType::assembly});

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
    .vtype = ipu::batchaffine::VertexType::cpp,
    .fillAlgo = ipu::batchaffine::partition::Algorithm::roundRobin
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
    .vtype = ipu::batchaffine::VertexType::multi,
    .fillAlgo = ipu::batchaffine::partition::Algorithm::roundRobin
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
    .vtype = ipu::batchaffine::VertexType::multiasm,
    .fillAlgo = ipu::batchaffine::partition::Algorithm::roundRobin
  });

  driver.compare_local(queries, refs);
  auto aln_results = driver.get_result();
  checkResults(aln_results);
}

TEST(PartitioningTest, FillFull) {
  int tilesUsed = 10;
  int maxBatches = 2;
  int maxAB = 10;
  int bufsize = maxBatches * maxAB;
  auto driver = ipu::batchaffine::SWAlgorithm({}, {
    .tilesUsed = tilesUsed,
    .maxAB = maxAB,
    .maxBatches = maxBatches,
    .bufsize = bufsize,
    .vtype = ipu::batchaffine::VertexType::cpp,
    .fillAlgo = ipu::batchaffine::partition::Algorithm::fillFirst
  });

  std::vector<std::string> a, b;
  for (int i = 0; i < tilesUsed * maxBatches; ++i) {
    a.push_back(std::string("A", maxAB));
    b.push_back(std::string("A", maxAB));
  }

  int errval = 0;
  auto buckets = driver.fillBuckets(a, b, errval);
  EXPECT_EQ(errval, 0);
}

TEST(PrepareTest, simple) {
  int tilesUsed = 10;
  int maxBatches = 2;
  int maxAB = 10;
  int bufsize = maxBatches * maxAB;
  ipu::batchaffine::IPUAlgoConfig config = {
    .tilesUsed = tilesUsed,
    .maxAB = maxAB,
    .maxBatches = maxBatches,
    .bufsize = bufsize,
    .vtype = ipu::batchaffine::VertexType::cpp,
    .fillAlgo = ipu::batchaffine::partition::Algorithm::fillFirst
  };

  std::vector<std::string> A, B;
  for (int i = 0; i < tilesUsed * maxBatches; ++i) {
    A.push_back(std::string("A", maxAB));
    B.push_back(std::string("A", maxAB));
  }

  std::vector<int32_t> inputs(config.getInputBufferSize32b() + 2);
  std::vector<int> mapping;

  inputs[0] = 0xDEADBEEF;
  inputs[inputs.size() - 1] = 0xDEADBEEF;
  ipu::batchaffine::SWAlgorithm::prepare_remote(config, A, B, &*inputs.begin() + 1, &*inputs.end() - 1, mapping);
  std::vector<int32_t> slice_begin(inputs.begin(), inputs.begin() + 10);
  std::vector<int32_t> slice_end(inputs.end() - 10, inputs.end());
  EXPECT_EQ(*inputs.begin(), 0xDEADBEEF) << "Start overwritten: " << swatlib::printVector(slice_begin);
  EXPECT_EQ(*(inputs.end() - 1), 0xDEADBEEF) << "End overwritten: " << swatlib::printVector(slice_end);
}