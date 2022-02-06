#include <plog/Log.h>
#undef LOG

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>

#include "ssw.hpp"
#include "gtest/gtest.h"

#include "ipu_batch_affine.h"
#include "swatlib/vector.hpp"

using std::max;
using std::min;
using std::string;
using std::vector;
using std::tuple;

using namespace StripedSmithWaterman;
vector<tuple<string, string>> INPUT_BATCHS = {
  {"/global/D1/projects/ipumer/inputs_ab/batch_0_A.txt", "/global/D1/projects/ipumer/inputs_ab/batch_0_B.txt"},
  {"/global/D1/projects/ipumer/inputs_ab/batch_1_A.txt", "/global/D1/projects/ipumer/inputs_ab/batch_1_B.txt"},
  {"/global/D1/projects/ipumer/inputs_ab/batch_2_A.txt", "/global/D1/projects/ipumer/inputs_ab/batch_2_B.txt"},
};

string aln2string(Alignment &aln) {
  std::stringstream ss;
  ss << "score=" << aln.sw_score << " score2=" << aln.sw_score_next_best;
  ss << " rbegin=" << aln.ref_begin << " rend=" << aln.ref_end;
  ss << " qbegin=" << aln.query_begin << " qend=" << aln.query_end;
  ss << " rend2=" << aln.ref_end_next_best << " mismatches=" << aln.mismatches;
  ss << " cigarstr=" << aln.cigar_string;
  return ss.str();
}

std::vector<std::string> loadSequences(const std::string& path) {
  std::vector<std::string> sequences;
  std::ifstream seqFile(path);
  string line;
  while (std::getline(seqFile, line)) {
    sequences.push_back(line);
  }
  return sequences;
}

class PerformanceBase : public ::testing::Test {
protected:
  vector<string> refs, queries;
};

class AlgoPerformance : public PerformanceBase, public ::testing::WithParamInterface<ipu::batchaffine::VertexType> {
};

TEST_P(AlgoPerformance, StressTest) {
  auto algotype = GetParam();
  int numWorkers = 8832;
  int numCmps = 40;
  int strlen = 150;
  if (algotype == ipu::batchaffine::VertexType::multi || algotype == ipu::batchaffine::VertexType::multiasm) {
    numWorkers = numWorkers / 6;
    numCmps = numCmps * 6;
  }

  // generate input strings
  int invalidRuns = 0;
  int numBatches = 10;
  auto driver = ipu::batchaffine::SWAlgorithm({}, {numWorkers, strlen, numCmps, numCmps * strlen * 2, algotype});
  for (int n = 0; n < numBatches; ++n) {
    refs = {};
    queries = {};
    for (int i = 0; i < numCmps * numWorkers; ++i) {
      refs.push_back(string(strlen, 'A'));
      queries.push_back(string(strlen, 'T'));
    }
    driver.compare_local(queries, refs, false);
    auto results = driver.get_result();
    bool valid = true;
    for (int i = 0; i < numCmps * numWorkers; ++i) {
      if (results.scores[i] != 0) {
        valid = false;
        // std::cout << driver.printTensors() << "\n";
        // exit(1);
      }
      // EXPECT_EQ(results.scores[i], 0) << "n: " << n << " mismatching score";
    }
    if (!valid) {
      PLOGW << "Invalid run encountered";
      invalidRuns++;
    }
  }
  EXPECT_EQ(invalidRuns, 0) << " invalid runs encountered in total " << numBatches << " runs";
}

TEST_P(AlgoPerformance, RunOptimalWithReverse) {
  auto algotype = GetParam();
  int numWorkers = 8832;
  int numCmps = 40;
  int strlen = 150;
  if (algotype == ipu::batchaffine::VertexType::multi || algotype == ipu::batchaffine::VertexType::multiasm) {
    numWorkers = numWorkers / 6;
    numCmps = numCmps * 6;
  }
  auto driver = ipu::batchaffine::SWAlgorithm({}, {numWorkers, strlen, numCmps, numCmps * strlen * 2, algotype});

  // generate input strings
  int numBatches = 10;
  for (int n = 0; n < numBatches; ++n) {
    refs = {};
    queries = {};
    for (int i = 0; i < numCmps * numWorkers; ++i) {
      refs.push_back(string(strlen, 'A'));
      queries.push_back(string(strlen, 'A'));
    }
    driver.compare_local(queries, refs, false);
    auto results = driver.get_result();
    for (int i = 0; i < numCmps * numWorkers; ++i) {
      EXPECT_EQ(results.scores[i], strlen) << "n: " << n << " mismatching score";
    }
  }
  // auto alns = driver.get_result();
  // for (int i = 0; i < numWorkers * numCmps; ++i) {
  //   std::cout << alns.scores[i] << " ";
  // }
  // std::cout << "\n";
}

TEST_P(AlgoPerformance, RunOptimal) {
  auto algotype = GetParam();
  int numWorkers = 8832;
  int numCmps = 40;
  int strlen = 150;
  if (algotype == ipu::batchaffine::VertexType::multi || algotype == ipu::batchaffine::VertexType::multiasm) {
    numWorkers = numWorkers / 6;
    numCmps = numCmps * 6;
  }
  auto driver = ipu::batchaffine::SWAlgorithm({}, {numWorkers, strlen, numCmps, numCmps * strlen * 2, algotype});

  // generate input strings
  int numBatches = 10;
  for (int n = 0; n < numBatches; ++n) {
    refs = {};
    queries = {};
    for (int i = 0; i < numCmps * numWorkers; ++i) {
      refs.push_back(string(strlen, 'A'));
      queries.push_back(string(strlen, 'T'));
    }
    driver.compare_local(queries, refs, false);
    auto results = driver.get_result();
    for (int i = 0; i < numCmps * numWorkers; ++i) {
      EXPECT_EQ(results.scores[i], 0) << "n: " << n << " mismatching score";
    }
  }
  // auto alns = driver.get_result();
  // for (int i = 0; i < numWorkers * numCmps; ++i) {
  //   std::cout << alns.scores[i] << " ";
  // }
  // std::cout << "\n";
}

INSTANTIATE_TEST_SUITE_P(
  VertexTypePerformance,
  AlgoPerformance,
  testing::Values(ipu::batchaffine::VertexType::cpp, ipu::batchaffine::VertexType::assembly, ipu::batchaffine::VertexType::multi, ipu::batchaffine::VertexType::multiasm)
  );

class PartitionPerformance : public PerformanceBase, public ::testing::WithParamInterface<ipu::partition::Algorithm> {
};

TEST_P(PartitionPerformance, RealBatches) {
  int numWorkers = 8832 / 6;
  int numCmps = 100;
  int strlen = 200;

  auto driver = ipu::batchaffine::SWAlgorithm({}, {numWorkers, strlen, numCmps, 8000 * 2, ipu::batchaffine::VertexType::multiasm, GetParam()});
  for (auto& [path_a, path_b] : INPUT_BATCHS) {
    refs = loadSequences(path_a);
    queries = loadSequences(path_b);
    // std::cout << "Len A: " << refs.size() << " Len B: " << queries.size() << "\n";
    driver.compare_local(refs, queries);
  }
}

INSTANTIATE_TEST_SUITE_P(
  PartitionTests,
  PartitionPerformance,
  testing::Values(ipu::partition::Algorithm::fillFirst, ipu::partition::Algorithm::roundRobin, ipu::partition::Algorithm::greedy)
  );