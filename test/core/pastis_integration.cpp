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
  auto batches = algo.create_batches({seqH, seqV}, {{0, 1, 1, 1}});
  auto job = algo.async_submit(&batches[0]);
  algo.blocking_join(*job);
  return job->batch->get_result().scores[0];
}

class PastisTest : public ::testing::Test {
protected:
  std::vector<ipu::Comparison> comps;
  std::vector<std::string> seqs;
  void SetUp() override {
    std::vector<std::string> seqH;
    loadSequences("/home/lukb/git/PASTIS-xdrop/build/seqs_seed_H.txt", seqH);
    std::vector<std::string> seqV;
    loadSequences("/home/lukb/git/PASTIS-xdrop/build/seqs_seed_V.txt", seqV);

    std::vector<size_t> seedseqH;
    loadSeedpos("/home/lukb/git/PASTIS-xdrop/build/seeds_H.txt", seedseqH);
    std::vector<size_t> seedseqV;
    loadSeedpos("/home/lukb/git/PASTIS-xdrop/build/seeds_V.txt", seedseqV);
    ASSERT_EQ(seqH.size(), seqV.size());
    ASSERT_EQ(seedseqH.size(), seedseqV.size());
    ASSERT_EQ(seqH.size(), seedseqH.size());

    seqs.resize(seqH.size() * 2);
    comps.resize(seqH.size());
    for (size_t i = 0; i < seqH.size(); i++) {
      seqs[i*2] = seqV[i];
      seqs[i*2+1] = seqH[i];
      comps[i] = {
          (int32_t)(2 * i), (int32_t)(2 * i + 1),
          seedseqV[i], seedseqH[i]
      };
    }
  }
};

TEST_F(PastisTest, XDropToSSW) {
    const ipu::SWConfig SW_CONFIGURATION = {
        .gapInit = 1,
        .gapExtend = -1,
        .matchValue = 1,
        .mismatchValue = -1,
        .ambiguityValue = -1,
        .similarity = swatlib::Similarity::blosum62,
        .datatype = swatlib::DataType::aminoAcid,
    };

    const ipu::IPUAlgoConfig ALGOCONFIG = {
        .numVertices = 1472,
        .maxSequenceLength = 1028,
        .maxComparisonsPerVertex = 360,
        .vertexBufferSize = 160000,
        .vtype = ipu::VertexType::xdropseedextend,
        .fillAlgo = ipu::Algorithm::fillFirst,
        .forwardOnly = false,
        .ioTiles = 0,
        .xDrop = 10,
        .bandPercentageXDrop = 0.5,
        .seedLength = 6,
    };


    auto encoder = getEncoder(swatlib::DataType::aminoAcid);
    auto sim = selectMatrix(swatlib::Similarity::blosum62, 1, -1, -1);
    int x = 0;
    for (size_t i = 0; i < comps.size(); i++) {
      auto encoded_seqH = encoder.encode(seqs[comps[i].indexA]);
      auto encoded_seqV = encoder.encode(seqs[comps[i].indexB]);
      PLOGE << seqs[comps[i].indexA]; 
      PLOGE << seqs[comps[i].indexB]; 
      PLOGE << seqs[comps[i].indexA].size() << "," << comps[i].seedAStartPos; 
      PLOGE << seqs[comps[i].indexB].size() << "," << comps[i].seedBStartPos; 
      auto score = ipumacore::xdrop::seed_extend_cpu<10, 1>(encoded_seqH, comps[i].seedAStartPos, encoded_seqV, comps[i].seedBStartPos, 6, sim);
      PLOGW << score;
      x += score;
    }
    
    auto algo = ipu::batchaffine::SWAlgorithm(SW_CONFIGURATION, ALGOCONFIG, 0, 1, 1);
    auto batches = algo.create_batches(seqs, comps);
    for (size_t i = 0; i < batches.size(); i++) {
      auto job = algo.async_submit(&batches[i]);
      PLOGW << "Wait for batch";
      algo.blocking_join(*job);
      // ASSERT_EQ(ipumacore::xdrop::seed_extend_cpu<>());
    }

    // for (size_t j = 0; j < batches.size(); j++) {
    //   auto &batch = batches[j];
    //   ipu::BlockAlignmentResults res = batch.get_result();
    //   // #pragma omp for
    //   for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
    //     auto orig_i = batch.origin_comparison_index[i];
    //     if (orig_i >= 0) {
    //       xscores[batch.origin_comparison_index[i]] = res.scores[i];
    //     }
    //   }
    // }
    
    ASSERT_EQ(10, x);
    // job->batch->get_result().scores[0];



  const int drop = 1000;
  // for(int i = 0; i <  seqHs.size(); i++) {
  //   auto max_len = std::max(seqHs[i].size(), seqVs[i].size()) + 1;
  //   std::string global_seqHs = initerstring(max_len) + seqHs[i];
  //   std::string global_seqVs = initerstring(max_len) + seqVs[i];

  //   auto ssw_score = alignFullSSW(global_seqHs, global_seqVs);

  //   auto encoder = getEncoder(DataType::nucleicAcid);
  //   auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
  //   auto ref = encoder.encode(global_seqVs);
  //   auto quer = encoder.encode(global_seqHs);
  //   auto core_score = alignSingleOnIpu(algo, global_seqHs, global_seqVs);

  //   // ASSERT_EQ(ssw_score, core_score) << "\nHSeq: " << seqHs[i] << "\nVSeq: " << seqVs[i] << alignSeqSeqan(seqHs[i], seqVs[i], drop);
  //   // ASSERT_GT(ssw_score/(double)core_score, 0.85);
  // }
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