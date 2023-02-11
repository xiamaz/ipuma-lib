#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <cmath>

#include "gtest/gtest.h"
#include "ssw/ssw.hpp"

// #define PRINT_DEBUG
#include "./helper.hpp"
#include "core/xdrop.hpp"

class ParityTest : public ::testing::Test {
 protected:
  void SetUp() override {
    seqHs = {
        "ATA"
        "AAA",
        "AAATA",
        "AAATA",
        "AAAA",
        "ATATCATCACTCCGATGGACG",
        "ATATCATCACTCCGATGGACGTTTCGATTATCGGCTGCGTGGGGAATGGCCCAGGTGAGGCGCTGGTTTCTACACTCGGCGTCACCGGCGGGACGAAACA",
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
        "AAA"
        "AAA",
        "AAATA",
        "AAAA",
        "AA",
        "TGGTTTCTACACTCGGCGT",
        "TGGTTTCTACACTCGGCGTCACCGGCGGCAACAAGAA",
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

    seedExtendSeqHs = {
        {2, 1, "ABBA"},
        {2, 0, "BBA"},
        {2, 1, "ABB"},
        {2, 1, "ABBA"}};
    seedExtendSeqVs = {
        {2, 0, "BBA"},
        {2, 1, "ABBA"},
        {2, 1, "ABB"},
        {2, 1, "ABBA"}};
    ASSERT_EQ(seqHs.size(), seqVs.size());
  }
  std::vector<std::string> seqHs, seqVs;
  std::vector<std::tuple<int, int, std::string>> seedExtendSeqHs, seedExtendSeqVs;
};

TEST_F(ParityTest, SeqAnToSSW) {
  for (int i = 0; i < seqHs.size(); i++) {
    auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
    std::string global_seqHs = initerstring(max_len) + seqHs[i] + initerstring(max_len);
    std::string global_seqVs = initerstring(max_len) + seqVs[i] + initerstring(max_len);

    auto seqan_score = alignSeqSeqan(global_seqHs, global_seqVs, 10000);
    auto ssw_score = alignFullSSW(global_seqHs, global_seqVs);
    ASSERT_EQ(seqan_score, ssw_score);
  }
}

TEST_F(ParityTest, XDropToSSW) {
  const int drop = 1000;
  for (int i = 0; i < seqHs.size(); i++) {
    auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
    std::string global_seqHs = initerstring(max_len) + seqHs[i];
    std::string global_seqVs = initerstring(max_len) + seqVs[i];

    auto ssw_score = alignFullSSW(global_seqHs, global_seqVs);

    auto encoder = getEncoder(DataType::nucleicAcid);
    auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
    auto ref = encoder.encode(global_seqVs);
    auto quer = encoder.encode(global_seqHs);
    auto core_score = ipumacore::xdrop::xdrop_doubleband_cpu<drop, 1, false>(ref, quer, sim);
    ASSERT_EQ(ssw_score, core_score) << "\nHSeq: " << seqHs[i] << "\nVSeq: " << seqVs[i] << alignSeqSeqan(seqHs[i], seqVs[i], drop);
    // ASSERT_GT(ssw_score/(double)core_score, 0.85);
  }
}

TEST_F(ParityTest, XDropToSeqAn) {
  const int drop = 10000;
  for (int i = 0; i < seqHs.size(); i++) {
    auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
    std::string global_seqHs = seqHs[i];
    std::string global_seqVs = seqVs[i];

    auto encoder = getEncoder(DataType::nucleicAcid);
    auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
    auto encoded_seqH = encoder.encode(global_seqHs);
    auto encoded_seqV = encoder.encode(global_seqVs);

    auto seqan_score = alignSeqSeqan(global_seqHs, global_seqVs, drop);
    auto core_score = ipumacore::xdrop::xdrop_doubleband_cpu<drop, 1, false>(encoded_seqH, encoded_seqV, sim);
    ASSERT_EQ(seqan_score, core_score);
  }
}

TEST_F(ParityTest, XDropToXDropRestrictedUnrestrictedX2) {
  const int drop = 2;
  for (int i = 0; i < seqHs.size(); i++) {
    auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
    std::string global_seqHs = seqHs[i];
    std::string global_seqVs = seqVs[i];

    auto encoder = getEncoder(DataType::nucleicAcid);
    auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
    auto encoded_seqH = encoder.encode(global_seqHs);
    auto encoded_seqV = encoder.encode(global_seqVs);

    auto core_score = ipumacore::xdrop::xdrop_doubleband_cpu<drop, 1, false>(encoded_seqH, encoded_seqV, sim);
    auto core_restricted_score = ipumacore::xdrop::xdrop_doubleband_restricted_cpu<drop, 1, 1000>(encoded_seqH, encoded_seqV, sim);
    ASSERT_EQ(core_score, core_restricted_score) << "SeqH: " << global_seqHs << ", SeqV: " << global_seqVs;
  }
}

TEST_F(ParityTest, XDropToXDropRestrictedX2) {
  const int drop = 2;
  for (int i = 0; i < seqHs.size(); i++) {
    auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
    std::string global_seqHs = seqHs[i];
    std::string global_seqVs = seqVs[i];

    auto encoder = getEncoder(DataType::nucleicAcid);
    auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
    auto encoded_seqH = encoder.encode(global_seqHs);
    auto encoded_seqV = encoder.encode(global_seqVs);

    auto core_score = ipumacore::xdrop::xdrop_doubleband_cpu<drop, 1, false>(encoded_seqH, encoded_seqV, sim);
    auto core_restricted_score = ipumacore::xdrop::xdrop_doubleband_restricted_cpu<drop, 1, 1000>(encoded_seqH, encoded_seqV, sim);
    ASSERT_EQ(core_score, core_restricted_score) << "SeqH: " << global_seqHs << ", SeqV: " << global_seqVs;
  }
}

TEST_F(ParityTest, XDropSeedExtend) {
  const int drop = 20;
  for (int i = 0; i < seedExtendSeqVs.size(); i++) {
    ASSERT_EQ(std::get<0>(seedExtendSeqHs[i]), std::get<0>(seedExtendSeqVs[i]));
    std::string global_seqHs = std::get<2>(seedExtendSeqHs[i]);
    std::string global_seqVs = std::get<2>(seedExtendSeqVs[i]);

    auto encoder = getEncoder(DataType::nucleicAcid);
    auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
    auto encoded_seqH = encoder.encode(global_seqHs);
    auto encoded_seqV = encoder.encode(global_seqVs);

    auto extended_score = ipumacore::xdrop::seed_extend_cpu<drop, 1>(
        encoded_seqH, std::get<1>(seedExtendSeqHs[i]),
        encoded_seqV, std::get<1>(seedExtendSeqVs[i]),
        std::get<0>(seedExtendSeqVs[i]), sim);
    PLOGW << extended_score;
    // ASSERT_EQ(core_score, core_restricted_score) << "SeqH: " << global_seqHs << ", SeqV: " << global_seqVs;
  }
}