#include <cmath>
#include "gtest/gtest.h"

#include "ssw/ssw.hpp"
#include "ipu_batch_affine.h"

inline std::string aln2string(StripedSmithWaterman::Alignment &aln) {
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

  int matchScore = ALN_MATCH_SCORE;
  int mismatchScore = ALN_MISMATCH_COST;
  int gapOpening = ALN_GAP_OPENING_COST;
  int gapExtend = ALN_GAP_EXTENDING_COST;
  int ambiguityCost = ALN_AMBIGUITY_COST;

  std::vector<std::string> queries, refs;
  
  void checkResults(const std::vector<StripedSmithWaterman::Alignment>& alns_ipu) {

    StripedSmithWaterman::Aligner ssw_aligner(matchScore, mismatchScore, gapOpening, gapExtend, ambiguityCost);

    StripedSmithWaterman::Filter ssw_filter;
    std::vector<StripedSmithWaterman::Alignment> alns(queries.size());
    for (int i = 0; i < queries.size(); ++i) {
      auto reflen = refs[i].size();
      auto qlen = queries[i].size();

      ssw_aligner.Align(refs[i].data(), reflen, queries[i].data(), qlen, ssw_filter, &alns[i], std::max((int)(qlen / 2), 15));
      EXPECT_EQ(alns_ipu[i].sw_score, alns[i].sw_score) << i << ": IPU score result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].query_begin, alns[i].query_begin) << i << ": IPU reference start result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].query_end, alns[i].query_end) << i << ": IPU reference end result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].ref_begin, alns[i].ref_begin) << i << ": IPU query start result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].ref_end, alns[i].ref_end) << i << ": IPU query end result does not match CPU SSW";
    }
  }
};

void test_aligns_ipu(std::vector<StripedSmithWaterman::Alignment> &alns, std::vector<std::string> query, std::vector<std::string> ref, ipu::batchaffine::SWAlgorithm &algo) {
  alns.reserve(query.size());

  algo.compare_local(query, ref);
  auto aln_results = algo.get_result();

  for (int i = 0; i < query.size(); i++) {
    alns.push_back({});
    StripedSmithWaterman::Alignment &aln = alns[i];

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

void test_aligns_ipu_mn(std::vector<StripedSmithWaterman::Alignment> &alns, std::vector<std::string> query, std::vector<std::string> ref, ipu::batchaffine::SWAlgorithm &algo) {
  alns.reserve(query.size());

  ipu::RawSequences seqs;
  ipu::Comparisons cmps;
  for (int i = 0; i < query.size(); ++i) {
    seqs.push_back(query[i]);
    seqs.push_back(ref[i]);
    cmps.push_back({i*2, i*2+1});
  }

  algo.compare_mn_local(seqs, cmps);
  auto aln_results = algo.get_result();

  for (int i = 0; i < query.size(); i++) {
    alns.push_back({});
    StripedSmithWaterman::Alignment &aln = alns[i];

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

TEST_F(ParityTest, UseAssembly) {
  int numWorkers = 10;
  int numCmps = 10;
  int strlen = 120;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -(gapOpening-gapExtend),
    .gapExtend = -gapExtend,
    .matchValue = matchScore,
    .mismatchValue = -mismatchScore,
    .ambiguityValue = -ambiguityCost,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::assembly, ipu::Algorithm::greedy});

  std::vector<StripedSmithWaterman::Alignment> alns_ipu(queries.size());
  test_aligns_ipu(alns_ipu, refs, queries, driver);

  checkResults(alns_ipu);
}

TEST_F(ParityTest, UseMultiAssembly) {
  int numWorkers = 10;
  int numCmps = 10;
  int strlen = 120;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -(gapOpening-gapExtend),
    .gapExtend = -gapExtend,
    .matchValue = matchScore,
    .mismatchValue = -mismatchScore,
    .ambiguityValue = -ambiguityCost,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::multiasm, ipu::Algorithm::greedy});

  std::vector<StripedSmithWaterman::Alignment> alns_ipu(queries.size());
  test_aligns_ipu(alns_ipu, refs, queries, driver);

  checkResults(alns_ipu);
}

TEST_F(ParityTest, UseMNComparison) {
  int numWorkers = 10;
  int numCmps = 10;
  int strlen = 120;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -(gapOpening-gapExtend),
    .gapExtend = -gapExtend,
    .matchValue = matchScore,
    .mismatchValue = -mismatchScore,
    .ambiguityValue = -ambiguityCost,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::multiasm, ipu::Algorithm::greedy});

  std::vector<StripedSmithWaterman::Alignment> alns_ipu(queries.size());
  test_aligns_ipu_mn(alns_ipu, refs, queries, driver);

  checkResults(alns_ipu);
}

TEST_F(ParityTest, UseGapOpeningCostMultiAsm) {
  gapOpening = 3;
  int numWorkers = 10;
  int numCmps = 10;
  int strlen = 120;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -(gapOpening-gapExtend),
    .gapExtend = -gapExtend,
    .matchValue = matchScore,
    .mismatchValue = -mismatchScore,
    .ambiguityValue = -ambiguityCost,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::multiasm, ipu::Algorithm::greedy});

  std::vector<StripedSmithWaterman::Alignment> alns_ipu(queries.size());
  test_aligns_ipu(alns_ipu, refs, queries, driver);

  checkResults(alns_ipu);
}

class AAParityTest : public ParityTest {
protected:
  void SetUp() override {
    queries = {
      "AAAAAA",
      "MMMMMM",
      "MATGGRRGAAAAP",
      "MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHL",
    };
    refs = {
      "AAAAAA",
      "MMMMMM",
      "MATGGRRGAAAAP",
      "MATGGRRGAAAAPLLYPGEVCPGMDIRNNLVAVAALLLGAAGHLTRLHELENCSVIEGHL",
    };
    gapOpening = 3;
    gapExtend = 1;
  }

  void checkResults(const std::vector<StripedSmithWaterman::Alignment>& alns_ipu) {
	  const std::vector<int8_t> mat50 = {
    //  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
      	5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,	// A
       -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,	// R
       -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,	// N
       -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,	// D
       -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,	// C
       -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,	// Q
       -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,	// E
      	0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,	// G
       -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,	// H
       -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,	// I
       -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,	// L
       -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,	// K
       -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,	// M
       -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,	// F
       -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,	// P
      	1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,	// S
      	0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, 	// T
       -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, 	// W
       -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, 	// Y
      	0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, 	// V
       -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, 	// B
       -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, 	// Z
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, 	// X
       -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1 	// *
	  };

	  /* This table is used to transform amino acid letters into numbers. */
	  const std::vector<int8_t> aa_table = {
	  	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	  	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	  	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	  	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	  	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	  	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
	  	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	  	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
	  };

    int matSize = static_cast<int>(std::sqrt(mat50.size()));

    StripedSmithWaterman::Aligner ssw_aligner(mat50.data(), matSize, aa_table.data(), aa_table.size());
    StripedSmithWaterman::Filter ssw_filter;
    std::vector<StripedSmithWaterman::Alignment> alns(queries.size());
    for (int i = 0; i < queries.size(); ++i) {
      auto reflen = refs[i].size();
      auto qlen = queries[i].size();

      ssw_aligner.Align(refs[i].data(), reflen, queries[i].data(), qlen, ssw_filter, &alns[i], std::max((int)(qlen / 2), 15));
      EXPECT_EQ(alns_ipu[i].sw_score, alns[i].sw_score) << i << ": IPU score result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].query_begin, alns[i].query_begin) << i << ": IPU reference start result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].query_end, alns[i].query_end) << i << ": IPU reference end result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].ref_begin, alns[i].ref_begin) << i << ": IPU query start result does not match CPU SSW";
      EXPECT_EQ(alns_ipu[i].ref_end, alns[i].ref_end) << i << ": IPU query end result does not match CPU SSW";
    }
  }
};

TEST_F(AAParityTest, UseCppVertex) {
  int numWorkers = 1;
  int numCmps = 10;
  int strlen = 120;
  int bufsize = 1000;
  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -(gapOpening-gapExtend),
    .gapExtend = -gapExtend,
    .matchValue = matchScore,
    .mismatchValue = -mismatchScore,
    .ambiguityValue = -ambiguityCost,
    .similarity = swatlib::Similarity::blosum50,
    .datatype = swatlib::DataType::aminoAcid,
  }, {numWorkers, strlen, numCmps, bufsize, ipu::VertexType::cpp, ipu::Algorithm::greedy});

  std::vector<StripedSmithWaterman::Alignment> alns_ipu(queries.size());
  test_aligns_ipu(alns_ipu, refs, queries, driver);

  checkResults(alns_ipu);
}