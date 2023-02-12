#include "gtest/gtest.h"
#include <plog/Log.h>

#include "ipu_batch_affine.h"

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
  }
  ASSERT_EQ(validCmps, Cmps.size());
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