#include "gtest/gtest.h"

#include "ipu_batch_affine.h"

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