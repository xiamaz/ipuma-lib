#include "gtest/gtest.h"
#include <plog/Log.h>

#include "ipu_batch_affine.h"

TEST(BatchingTest, testOriginComparisonIndex) {
  int sequence_length = 100;
  int batches_count = 2;
  ipu::SWConfig swconfig;
  ipu::IPUAlgoConfig config;
  ipu::RawSequences Seqs;
  ipu::Comparisons Cmps;
  swconfig.datatype = swatlib::DataType::aminoAcid;
  config = {
    2 /*Tiles used*/,
    10000 /*maxSequenceLength*/,
    20,
    sequence_length * 2 * 10,
    ipu::VertexType::xdroprestrictedseedextend,
    ipu::Algorithm::greedy,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  auto totalCmpCount = numberCmps * batches_count;
  for (int i = 0; i < totalCmpCount; ++i) {
    Seqs.push_back(std::string(i + 1, 'A'));
    Seqs.push_back(std::string(i + 100, 'B'));
    Cmps.push_back({.originalComparisonIndex = i, .indexA = i * 2, .sizeA = (int32_t) Seqs[i*2].size(), .indexB = i * 2 + 1, .sizeB = (int32_t) Seqs[i*2+1].size(), .seeds = {{{10, 10}, {100, 100}}}, .complexity = 0});
  }

  auto batches = ipu::create_batches(Seqs, Cmps, config, swconfig);

  int validCmps = 0;
  uint64_t cellCount = 0;
  std::map<size_t, int> hasCmp;
  for (auto batch : batches) {
    ipu::XDropMeta* m = reinterpret_cast<ipu::XDropMeta*>(batch.getMetaBuffer());
    for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
      auto orig_i_packed = batch.origin_comparison_index[i];
      auto [orig_i, seed_j] = ipu::unpackOriginIndex(orig_i_packed);
      if (orig_i >= 0) {
        auto bm = m[i];
        auto sizeA = bm.sizeA;
        auto sizeB = bm.sizeB;
        auto iA = sizeA - 1;
        auto iB = sizeB - 100;
        ASSERT_EQ(iA, orig_i);
        ASSERT_EQ(iA, iB);
        if (hasCmp.find(iA) != hasCmp.end()) {
          hasCmp[iA] += 1;
        } else {
          hasCmp[iA] = 1;
        }
      }
    }
    cellCount += batch.cellCount;
  }
  for (const auto& [k, v] : hasCmp) {
    ASSERT_EQ(v, NSEEDS);
    ASSERT_LE(k, totalCmpCount);
  }
  for (int i = 0; i < totalCmpCount; ++i) {
    ASSERT_NE(hasCmp.find(i), hasCmp.end());
    ASSERT_EQ(hasCmp[i], NSEEDS);
  }
}

TEST(BatchingTest, testValidNumberCmps) {
  int sequence_length = 100;
  int batches_count = 2;
  ipu::SWConfig swconfig;
  ipu::IPUAlgoConfig config;
  ipu::RawSequences Seqs;
  ipu::Comparisons Cmps;
  swconfig.datatype = swatlib::DataType::aminoAcid;
  config = {
    2 /*Tiles used*/,
    100 /*maxSequenceLength*/,
    20,
    sequence_length * 2 * 10,
    ipu::VertexType::xdroprestrictedseedextend,
    ipu::Algorithm::fillFirst,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  for (int i = 0; i < numberCmps * batches_count; ++i) {
    Seqs.push_back(std::string(sequence_length, 'A'));
    Seqs.push_back(std::string(sequence_length, 'B'));
    Cmps.push_back({.originalComparisonIndex = i, .indexA = i * 2, .sizeA = (int32_t) Seqs[i*2].size(), .indexB = i * 2 + 1, .sizeB = (int32_t) Seqs[i*2+1].size(), .seeds = {{{10, 10}, {100, 100}}}, .complexity = 0});
  }

  auto batches = ipu::create_batches(Seqs, Cmps, config, swconfig);

  int validCmps = 0;
  uint64_t cellCount = 0;
  for (auto batch : batches) {
    ipu::XDropMeta* m = reinterpret_cast<ipu::XDropMeta*>(batch.getMetaBuffer());
    for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
      auto orig_i_packed = batch.origin_comparison_index[i];
      auto [orig_i, seed_j] = ipu::unpackOriginIndex(orig_i_packed);
      if (orig_i >= 0) {
        auto orig_cmp = Cmps[orig_i];
        auto batch_meta = m[i];
        validCmps += 1;
      }
    }
    cellCount += batch.cellCount;
  }
  ASSERT_EQ(validCmps, Cmps.size() * NSEEDS);
  ASSERT_EQ(cellCount, numberCmps * batches_count * sequence_length * sequence_length * NSEEDS);
}

TEST(BatchingTest, testUniformCorrectness) {
  int sequence_length = 100;
  int batches_count = 1;
  ipu::SWConfig swconfig;
  ipu::IPUAlgoConfig config;
  ipu::RawSequences Seqs;
  ipu::Comparisons Cmps;
  swconfig.datatype = swatlib::DataType::aminoAcid;
  config = {
    1 /*Tiles used*/,
    100 /*maxSequenceLength*/,
    14,
    sequence_length * 2 * 10,
    ipu::VertexType::xdroprestrictedseedextend,
    ipu::Algorithm::fillFirst,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  for (int i = 0; i < numberCmps * batches_count; ++i) {
    Seqs.push_back(std::string(sequence_length, 'A'));
    Seqs.push_back(std::string(sequence_length, 'B'));
    Cmps.push_back({.originalComparisonIndex = i, .indexA = i * 2, .sizeA = (int32_t) Seqs[i*2].size(), .indexB = i * 2 + 1, .sizeB = (int32_t) Seqs[i*2+1].size(), .seeds = {{{10, 10}, {100, 100}}}, .complexity = 0});
  }

  auto batches = ipu::create_batches(Seqs, Cmps, config, swconfig);

  for (auto batch : batches) {
    PLOGI << batch.toString();
    ipu::XDropMeta* m = reinterpret_cast<ipu::XDropMeta*>(batch.getMetaBuffer());
    for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
      auto orig_i_packed = batch.origin_comparison_index[i];
      auto [orig_i, seed_j] = ipu::unpackOriginIndex(orig_i_packed);
      if (orig_i >= 0) {
        auto orig_cmp = Cmps[orig_i];
        auto batch_meta = m[i];
        ASSERT_EQ(Seqs[orig_cmp.indexA].size(), batch_meta.sizeA);
        ASSERT_EQ(Seqs[orig_cmp.indexB].size(), batch_meta.sizeB);
        ASSERT_EQ(orig_cmp.seeds[seed_j].seedAStartPos, batch_meta.seedAStartPos);
        ASSERT_EQ(orig_cmp.seeds[seed_j].seedBStartPos, batch_meta.seedBStartPos);
      }
    }
  }
}