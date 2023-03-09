#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <omp.h>

#include <cmath>
#include <iostream>
#include <string_view>

#include "gtest/gtest.h"
#include "ssw/ssw.hpp"

// #define PRINT_DEBUG
#include "./helper.hpp"
#include "core/xdrop.hpp"
#include "ipu_batch_affine.h"

// int alignSingleOnIpu(ipu::batchaffine::SWAlgorithm& algo, const std::string& seqH, const std::string& seqV) {
//   auto batches = algo.create_batches({seqH, seqV}, {{0, 1, {{1, 1}}}});
//   auto job = algo.async_submit(&batches[0]);
//   algo.blocking_join(*job);
//   return job->batch->get_result().scores[0][0];
// }

const std::string basedir = "/home/lukb/git/ELBA/ELBA/build_release/";
// const std::string basedir = "/global/D1/projects/ipumer/inputs_ab/";

class ElbaTest : public ::testing::Test {
 protected:
  std::vector<ipu::Comparison> comps;
  std::vector<std::string_view> seqs;
  std::vector<std::string> seqH;
  std::vector<std::string> seqV;
  void SetUp() override {
    loadSequences(basedir + "seqs_seed_H.txt", seqH);
    loadSequences(basedir + "seqs_seed_V.txt", seqV);

    std::vector<size_t> seedseqH1;
    loadSeedpos(basedir + "seeds_H1.txt", seedseqH1);
    std::vector<size_t> seedseqV1;
    loadSeedpos(basedir + "seeds_V1.txt", seedseqV1);

    std::vector<size_t> seedseqH2;
    loadSeedpos(basedir + "seeds_H2.txt", seedseqH2);
    std::vector<size_t> seedseqV2;
    loadSeedpos(basedir + "seeds_V2.txt", seedseqV2);

    ASSERT_EQ(seqH.size(), seqV.size());
    ASSERT_EQ(seedseqH1.size(), seedseqV1.size());
    ASSERT_EQ(seedseqH2.size(), seedseqV2.size());
    ASSERT_EQ(seedseqH1.size(), seedseqH2.size());
    ASSERT_EQ(seqH.size(), seedseqH1.size());

    seqs.resize(seqH.size() * 2);
    // comps.resize(seqH.size());
    comps.resize(seqH.size());

    std::vector<size_t> complexities;
    loadSeedpos("complexity.txt", complexities);
    for (int i{}; i< 10; i++) {
      PLOGD << complexities[i];
    }

    // comps.resize(200);
    for (size_t i = 0; i < comps.size(); i++) {
      // if (i > 9660 && i < 10000) {
      //   std::string A = seqV[i].substr(seedseqV[i], 17);
      //   std::string B = seqH[i].substr(seedseqH[i], 17);
      //   PLOGW << "=======" << i << "========";
      //   PLOGW << A << ": " << seqV[i].size() << "," << seedseqV[i];
      //   PLOGW << B << ": " << seqH[i].size() << "," << seedseqH[i];
      //   PLOGW << std::endl;
      // }
      // PLOGW.printf("%d::%d, %d::%d", (int32_t)seedseqV1[i], (int32_t)seedseqH1[i], (int32_t)seedseqV2[i], (int32_t)seedseqH2[i]);
      seqs[i * 2] = seqV[i];
      seqs[i * 2 + 1] = seqH[i];
      comps[i] = {
          (int64_t) i,
          (int32_t)(2 * i),
          seqs[i * 2].size(),
          (int32_t)(2 * i + 1),
          seqs[i * 2 + 1].size(),
          {{{(int32_t)seedseqV1[i], (int32_t)seedseqH1[i]},
            {(int32_t)seedseqV2[i], (int32_t)seedseqH2[i]}}},
          complexities[i],
      };
    }
  }
};

TEST_F(ElbaTest, XDropRun) {
  const int X = 10;
  // const ipu::SWConfig SW_CONFIGURATION = {
  //     .gapInit = 1,
  //     .gapExtend = -1,
  //     .matchValue = 1,
  //     .mismatchValue = -1,
  //     .ambiguityValue = -1,
  //     .similarity = swatlib::Similarity::blosum62,
  //     .datatype = swatlib::DataType::aminoAcid,
  // };

  // const ipu::IPUAlgoConfig ALGOCONFIG = {
  //     .numVertices = 1472,
  //     .maxSequenceLength = 1028,
  //     .maxComparisonsPerVertex = 360,
  //     .vertexBufferSize = 160000,
  //     .vtype = ipu::VertexType::xdropseedextend,
  //     .fillAlgo = ipu::Algorithm::fillFirst,
  //     .forwardOnly = false,
  //     .ioTiles = 0,
  //     .xDrop = 10,
  //     .bandPercentageXDrop = 0.5,
  //     .seedLength = 6,
  // };

  const ipu::SWConfig SW_CONFIGURATION = {
      .gapInit = -1,
      .gapExtend = -1,
      .matchValue = 1,
      .mismatchValue = -1,
      .ambiguityValue = -1,
      .similarity = swatlib::Similarity::nucleicAcid,
      .datatype = swatlib::DataType::nucleicAcid,
  };

  const ipu::IPUAlgoConfig ALGOCONFIG = {
      .numVertices = 1472,
      .maxSequenceLength = 19295,
      .maxComparisonsPerVertex = 10 * 2,
      .vertexBufferSize = 160000,
      .vtype = ipu::VertexType::xdroprestrictedseedextend,
      .fillAlgo = ipu::Algorithm::greedy,
      .forwardOnly = false,
      .ioTiles = 0,
      .xDrop = X,
      .bandPercentageXDrop = 0.45,
      .seedLength = 17,
  };

  auto encoder = getEncoder(swatlib::DataType::nucleicAcid);
  auto sim = selectMatrix(swatlib::Similarity::nucleicAcid, 1, -1, -1);
  int x = 0;
  // for (size_t i = 0; i < comps.size(); i++) {
  // {
  //   size_t i = 627;
  //   std::cout << i << std::endl;
  //   auto encoded_seqH = encoder.encode(seqs[comps[i].indexA]);
  //   auto encoded_seqV = encoder.encode(seqs[comps[i].indexB]);
  //   PLOGE << seqs[comps[i].indexA];
  //   PLOGE << seqs[comps[i].indexB];
  //   PLOGE << seqs[comps[i].indexA].size() << "," << comps[i].seedAStartPos;
  //   PLOGE << seqs[comps[i].indexB].size() << "," << comps[i].seedBStartPos;
  //   auto Asize = comps[i].seedAStartPos;
  //   auto Ain = encoded_seqH;

  //   auto Bsize = comps[i].seedBStartPos;
  //   auto Bin = encoded_seqV;

  //   int larger = Asize > Bsize;
  //   #define SWAP(x, a, b) ((x) ? (a) : (b))
  //   #define RSWAP(x, a, b) ((x) ? (a) : (b))

  //   auto score = ipumacore::xdrop::seed_extend_restricted_cpu<10, 1>(
  //     SWAP(larger, Ain, Bin),
  //     SWAP(larger, Asize, Bsize),
  //     RSWAP(larger, Ain, Bin),
  //     RSWAP(larger, Asize, Bsize),
  //     17,  sim, (int)(19295 * 0.45));
  //   // PLOGW << score;
  //   x += score;
  // }

  auto algo = ipu::batchaffine::SWAlgorithm(SW_CONFIGURATION, ALGOCONFIG, 0, 1, 1);
  std::cout << "START BATCHES" << std::endl;
  auto batches = algo.create_batches(seqs, comps);
  std::cout << "DONE BATCHES" << std::endl;


  std::vector<size_t> complexities;
  loadSeedpos("complexity.txt", complexities);



  std::vector<int> complexity(comps.size(), 0);

  for (size_t i = 0; i < batches.size(); i++) {
    std::cout << "SUBMIT" << std::endl;
    auto job = algo.async_submit(&batches[i]);
    PLOGW << "Wait for batch";
    algo.blocking_join(*job);
    // ASSERT_EQ(ipumacore::xdrop::seed_extend_cpu<>());
    int j = i;
    auto& batch = batches[j];
    ipu::BlockAlignmentResults res = batch.get_result();
    // #pragma omp parallel for
    for (int k = 0; k < batch.origin_comparison_index.size(); ++k) {
      auto orig_i = batch.origin_comparison_index[k];
      if (orig_i >= 0) {
        for (int seed = 0; seed < NSEEDS; seed++) {
          // PLOGW.printf("%d -> %d\n", k, batch.origin_comparison_index[k]);
          // printf("AAA: %d %d %d %d %d %d %d\n",
          //   seqs[comps[orig_i].indexA].size(),
          //   seqs[comps[orig_i].indexB].size(),
          //   comps[orig_i].seedAStartPos,
          //   comps[orig_i].seedBStartPos,
          //   res.scores[k],
          //   res.a_range_result[k],
          //   res.b_range_result[k]);
          // complexity[orig_i] = res.a_range_result[k / NSEEDS][seed] + res.b_range_result[k / NSEEDS][seed];
          // if (abs((double) (complexities[orig_i] - complexity[orig_i]))/ (double) complexities[orig_i] * 100.0 > 2.0) {
          //   printf("complexity mismatch %d %d\n", complexities[orig_i],complexity[orig_i]);
          // }
          // continue;

          std::string sA{seqs[comps[orig_i].indexA]};
          auto encoded_seqH = encoder.encode(sA);

          std::string sB{seqs[comps[orig_i].indexB]};
          auto encoded_seqV = encoder.encode(sB);
          auto seedAPos = comps[orig_i].seeds[seed].seedAStartPos;
          auto Ain = encoded_seqH;

          auto seedBPos = comps[orig_i].seeds[seed].seedBStartPos;
          auto Bin = encoded_seqV;

          if (seedAPos == -1) {continue;}

          const int ipu_score_combined = res.a_range_result[k / NSEEDS][seed] + res.b_range_result[k / NSEEDS][seed];
          const int ipu_score_simple = res.scores[k / NSEEDS][seed];
          const int ipu_score = ipu_score_combined;
          complexity[orig_i] += ipu_score;

          // auto score = ipumacore::xdrop::seed_extend_restricted_cpu<X, 1>(
          //     Ain,
          //     seedAPos,
          //     Bin,
          //     seedBPos,
          //     17, sim, (int)(19295 * 0.45));


              // if (ipu_score_combined != ipu_score_simple) {
              //     printf("Expected combined %d, got from IPU %d = %d + %d\n", res.scores[k / NSEEDS][seed], res.a_range_result[k / NSEEDS][seed] + res.b_range_result[k / NSEEDS][seed], res.a_range_result[k / NSEEDS][seed] , res.b_range_result[k / NSEEDS][seed]);
              // }
              // if (abs(ipu_score-score)/(double)score * 100.0 > 1.0) {
              //     #pragma omp critical
              //     {
              //       printf("Expected Score %d, got from IPU %d (%f %%)\n", score, ipu_score, abs(score-ipu_score)/ (double)score * 100.0);
              //       printf("Missing: %d => (%d|%d)\n", Ain.size(),  seedAPos, Ain.size() - 17 - seedAPos);
              //       printf("Missing: %d => (%d|%d)\n\n", Bin.size(),  seedBPos, Bin.size() - 17 - seedBPos);
              //     }
              // }
          // ASSERT_EQ(res.scores[k / NSEEDS][seed], score) << "lenA: " << sA.size() << ", lenB: " << sB.size();
        }
      }
    }
  }
  {
    using namespace std;
    ofstream myfile;
    myfile.open("complexity_reorder.txt");
    for (auto&& i : complexity) {
      myfile << i << '\n';
    }
    myfile.close();
  }

  // ASSERT_EQ(10, x);
  // job->batch->get_result().scores[0];

  // const int drop = 1000;
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