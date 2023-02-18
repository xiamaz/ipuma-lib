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
    ipu::VertexType::xdropseedextend,
    ipu::Algorithm::greedy,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  auto totalCmpCount = numberCmps * batches_count;
  for (int i = 0; i < totalCmpCount; ++i) {
    Seqs.push_back(std::string(i + 1, 'A'));
    Seqs.push_back(std::string(i + 100, 'B'));
    Cmps.push_back({.indexA = i * 2, .indexB = i * 2 + 1});
  }

  auto batches = ipu::create_batches(Seqs, Cmps, config, swconfig);

  int validCmps = 0;
  uint64_t cellCount = 0;
  std::map<size_t, bool> hasCmp;
  for (size_t i = 0; i < totalCmpCount; ++i) {
    hasCmp[i] = false;
  }
  for (auto batch : batches) {
    ipu::XDropMeta* m = reinterpret_cast<ipu::XDropMeta*>(batch.getMetaBuffer());
    for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
      auto orig_i = batch.origin_comparison_index[i];
      if (orig_i >= 0) {
        auto bm = m[i];
        auto sizeA = bm.sizeA;
        auto sizeB = bm.sizeB;
        auto iA = sizeA - 1;
        auto iB = sizeB - 100;
        ASSERT_EQ(iA, iB);
        hasCmp[iA] = true;
      }
    }
    cellCount += batch.cellCount;
  }
  for (const auto& [k, v] : hasCmp) {
    ASSERT_TRUE(v);
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
    ipu::VertexType::xdropseedextend,
    ipu::Algorithm::fillFirst,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  for (int i = 0; i < numberCmps * batches_count; ++i) {
    Seqs.push_back(std::string(sequence_length, 'A'));
    Seqs.push_back(std::string(sequence_length, 'B'));
    Cmps.push_back({.indexA = i * 2, .indexB = i * 2 + 1});
  }

  auto batches = ipu::create_batches(Seqs, Cmps, config, swconfig);

  int validCmps = 0;
  uint64_t cellCount = 0;
  for (auto batch : batches) {
    ipu::XDropMeta* m = reinterpret_cast<ipu::XDropMeta*>(batch.getMetaBuffer());
    for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
      auto orig_i = batch.origin_comparison_index[i];
      if (orig_i >= 0) {
        auto orig_cmp = Cmps[orig_i];
        auto batch_meta = m[i];
        validCmps += 1;
      }
    }
    cellCount += batch.cellCount;
  }
  ASSERT_EQ(validCmps, Cmps.size());
  ASSERT_EQ(cellCount, numberCmps * batches_count * sequence_length * sequence_length);
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
    ipu::VertexType::xdropseedextend,
    ipu::Algorithm::fillFirst,
  };
  size_t totalBufsize = config.vertexBufferSize * config.numVertices;
  size_t numberCmps = totalBufsize / (sequence_length * 2);
  for (int i = 0; i < numberCmps * batches_count; ++i) {
    Seqs.push_back(std::string(sequence_length, 'A'));
    Seqs.push_back(std::string(sequence_length, 'B'));
    Cmps.push_back({.indexA = i * 2, .indexB = i * 2 + 1});
  }

  auto batches = ipu::create_batches(Seqs, Cmps, config, swconfig);

  for (auto batch : batches) {
    PLOGI << batch.toString();
    ipu::XDropMeta* m = reinterpret_cast<ipu::XDropMeta*>(batch.getMetaBuffer());
    for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
      auto orig_i = batch.origin_comparison_index[i];
      if (orig_i >= 0) {
        auto orig_cmp = Cmps[orig_i];
        auto batch_meta = m[i];
        ASSERT_EQ(Seqs[orig_cmp.indexA].size(), batch_meta.sizeA);
        ASSERT_EQ(Seqs[orig_cmp.indexB].size(), batch_meta.sizeB);
        ASSERT_EQ(orig_cmp.seedAStartPos, batch_meta.seedAStartPos);
        ASSERT_EQ(orig_cmp.seedBStartPos, batch_meta.seedBStartPos);
      }
    }
  }
}