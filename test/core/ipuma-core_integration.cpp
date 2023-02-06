#include <cmath>
#include "gtest/gtest.h"

#include "ssw/ssw.hpp"

#include <seqan/align.h>
#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

// #define PRINT_DEBUG 
#include "core/xdrop.hpp"
#include "ipu_batch_affine.h"

#include "./helper.hpp"


int alignSingleOnIpu(ipu::batchaffine::SWAlgorithm& algo, const std::string& seqH, const std::string& seqV) {
  auto batches = algo.create_batches({seqH, seqV}, {{0, 1}});
  auto job = algo.async_submit(&batches[0]);
  algo.blocking_join(*job);
  return job->batch->get_result().scores[0];
}



class ParityIntegrationTest : public ::testing::Test {
protected:
  void SetUp() override {
  seqHs = {
    "ATATCATCACTCCGATGGACG",
    "ATATCATCACTCCGATGGACGTTTCGATTATCGGCTGCGTGGGGAATGGCCCAGGTGAGGCGCTGGTTTCTACACTCGGCGTCACCGGCGGGACGAAACA",
    "ATA"
    "AAA",
    "AAATA",
    "AAATA",
    "AAAA",
    "AAAATTTTCCCCGGGG",
    "AAAATTTTCCCCGGGG",
    "AAAATTTTCCCCGGGG",
    "AAAATTTTCCCCGGGG",
    "AAAATTTTCCCCGGGG",
    "AAAATTTTCCCCGGGG",
    "AAAATTTTCCCCGGGG",
    "AATGAGAATGATGTCGTTCGAAATTTGACCAGTCAAACCGCGGGCAATAAGGTCTTCGTTCAGGGCATAGACCTTAATGGGGGCATTACGCAGACTTTCA",
    "ATCTGGCAGGTAAAGATGAGCTCAACAAAGTGATCCAGCATTTTGGCAAAGGAGGCTTTGATGTGATTACTCGCGGTCAGGTGCCACCTAACCCGTCTGA",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGGTATTGGGAGCACG",
    "CCCCGCACCCGCAAGCCGCCGAGAAAAAAAGGATGAGGGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGAACA",
    "CGTGCTCCCAATACCCGCCAGGAATTAGCGCGCTTACGCTCAGAGCTGGCGGAAAAAGAGAGCCAGCAACTGGATGCGTATTTGCAGGCCTTGCGTAATC",
    "CCCCGCACCCGCAAGCCGCCGAGAAAAAAAGGATGAGGGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGAACA",
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
  seqVs = {
    "TGGTTTCTACACTCGGCGT",
    "TGGTTTCTACACTCGGCGTCACCGGCGGCAACAAGAA",
    "AAA"
    "AAA",
    "AAATA",
    "AAAA",
    "AA",
    "AAAATTTTCCCCGGGG",
    "AAAATTTTACCCGGGG",
    "AAAATTTTACCCCGGGG",
    "AAAATTTCCCCGGGG",
    "GCTAGCTAGCTAGCTA",
    "GCTAAAATTTTCCCCGGGG",
    "AAAATTTTCCCCGGGGACT",
    "AATGAGAATGATGTCNTTCNAAATTTGACCAGTCAAACCGCGGGCAATAAGGTCTTCGTTCAGGGCATAGACCTTAATGGGGGCATTACGCAGACTTTCA",
    "ATCTGGCAGGTAAAGATGAGCTCAACAAAGTGATCCAGCATTTTGGCAAAGGAGGCTTTGATGTGATTACTCGCGGTCAGGTGCCACCTAANNCGTCTGA",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGTTATTGGCAGACAG",
    "GCACCGTCCAGCCAACCGCCGAGAAGAAAAGAATGAGTGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGACAC",
    "CTGTCTGCCAATAACCGCCAGGAATTAGCGCGCTTACGCTCAGAGCTGGCGGAAAAAGAGAGCCAGCAACTGGATGCGTATTTGCAGGCCTTGCGTAATC",
    "GCACCGTCCAGCCAACCGCCGAGAAGAAAAGAATGAGTGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGACAC",
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
  ASSERT_EQ(seqHs.size(), seqVs.size());
  }
  std::vector<std::string> seqHs, seqVs;
};

TEST_F(ParityIntegrationTest, XDropToSSW) {
auto algo = ipu::batchaffine::SWAlgorithm({
    .gapInit = 0, .gapExtend = -1, .matchValue = 1, .mismatchValue = -1, .ambiguityValue = -1,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {1, 1000, 10, 10000, ipu::VertexType::multixdrop});


  const int drop = 1000;
  for(int i = 0; i <  seqHs.size(); i++) {
    auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
    std::string global_seqHs = initerstring(max_len) + seqHs[i];
    std::string global_seqVs = initerstring(max_len) + seqVs[i];

    auto ssw_score = alignFullSSW(global_seqHs, global_seqVs);

    auto encoder = getEncoder(DataType::nucleicAcid);
    auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
    auto ref = encoder.encode(global_seqVs);
    auto quer = encoder.encode(global_seqHs);
    auto core_score = alignSingleOnIpu(algo, global_seqHs, global_seqVs);
    ASSERT_EQ(ssw_score, core_score) << "\nHSeq: " << seqHs[i] << "\nVSeq: " << seqVs[i] << alignSeqSeqan(seqHs[i], seqVs[i], drop);
    // ASSERT_GT(ssw_score/(double)core_score, 0.85);
  }
}

// TEST_F(ParityIntegrationTest, XDropToSeqAnX3) {
//   const int drop = 3;
//   for(int i = 0; i <  seqHs.size(); i++) {
//     auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
//     std::string global_seqHs = initerstring(max_len) + seqHs[i] + initerstring(max_len);
//     std::string global_seqVs = initerstring(max_len) + seqVs[i] + initerstring(max_len);

//     auto encoder = getEncoder(DataType::nucleicAcid);
//     auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
//     auto encoded_seqH = encoder.encode(global_seqHs);
//     auto encoded_seqV = encoder.encode(global_seqVs);


//     auto seqan_score = alignSeqSeqan(global_seqHs, global_seqVs, drop);
//     auto core_score = ipumacore::xdrop::xdrop_doubleband_cpu<drop, 1, false>(encoded_seqH, encoded_seqV, sim);
//     ASSERT_EQ(seqan_score, core_score);
//   }
// }


// TEST_F(ParityIntegrationTest, XDropToXDropRestrictedUnrestrictedX2) {
//   const int drop = 2;
//   for(int i = 0; i <  seqHs.size(); i++) {
//     auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
//     std::string global_seqHs = seqHs[i];
//     std::string global_seqVs = seqVs[i];

//     auto encoder = getEncoder(DataType::nucleicAcid);
//     auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
//     auto encoded_seqH = encoder.encode(global_seqHs);
//     auto encoded_seqV = encoder.encode(global_seqVs);


//     auto core_score = ipumacore::xdrop::xdrop_doubleband_cpu<drop, 1, false>(encoded_seqH, encoded_seqV, sim);
//     auto core_restricted_score = ipumacore::xdrop::xdrop_doubleband_restricted_cpu<drop, 1, 1000>(encoded_seqH, encoded_seqV, sim);
//     ASSERT_EQ(core_score, core_restricted_score) << "SeqH: " << global_seqHs << ", SeqV: " << global_seqVs;
//   }
// }

// TEST_F(ParityIntegrationTest, XDropToXDropRestrictedX2) {
//   const int drop = 2;
//   for(int i = 0; i <  seqHs.size(); i++) {
//     auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
//     std::string global_seqHs = seqHs[i];
//     std::string global_seqVs = seqVs[i];

//     auto encoder = getEncoder(DataType::nucleicAcid);
//     auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
//     auto encoded_seqH = encoder.encode(global_seqHs);
//     auto encoded_seqV = encoder.encode(global_seqVs);


//     auto core_score = ipumacore::xdrop::xdrop_doubleband_cpu<drop, 1, false>(encoded_seqH, encoded_seqV, sim);
//     auto core_restricted_score = ipumacore::xdrop::xdrop_doubleband_restricted_cpu<drop, 1, 1000>(encoded_seqH, encoded_seqV, sim);
//     ASSERT_EQ(core_score, core_restricted_score) << "SeqH: " << global_seqHs << ", SeqV: " << global_seqVs;
//   }
// }