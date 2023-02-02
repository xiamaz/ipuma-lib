#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>

#include <algorithm>
#include <cxxopts.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <span>
#include <string>
#include <vector>

#include "./data.h"
#include "./xdrop_cpu.h"

#include <utility>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include "../src/ipuma.h"
#include "../src/driver.hpp"


using std::max;
using std::min;
using std::string;
using std::vector;
using std::tuple;

// const std::string basepath = "/Users/lbb/git/ipuma-lib/xdrop/data/inputs_ab/";
const std::string basepath = "/global/D1/projects/ipumer/inputs_ab/";
vector<tuple<string, string>> INPUT_BATCHS = {
  {basepath + "batch_0_A.txt", basepath + "batch_0_B.txt"},
  {basepath + "batch_1_A.txt", basepath + "batch_1_B.txt"},
  {basepath + "batch_2_A.txt", basepath + "batch_2_B.txt"},
};

// vector<tuple<string, string>> RR_ERROR_BATCHS = {
//   {basepath + "rr_As.txt", basepath + "rr_Bs.txt"},
// };

std::vector<std::string> loadSequences(const std::string& path) {
  std::vector<std::string> sequences;
  std::ifstream seqFile(path);
  string line;
  while (std::getline(seqFile, line)) {
    sequences.push_back(line);
  }
  return sequences;
}

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::verbose, &consoleAppender);

  //    0 1 2 3 4 5 6 7 n
  // 0 [0,0,0,0,0,0,0,0,0]
  // 1 [0,1,1,0,0,1,0,1,1]
  // 2 [0,1,2,1,0,1,0,1,2]
  // 3 [0,0,1,3,2,1,0,0,1]
  // 4 [0,0,0,2,4,3,2,1,0]
  // m [0,1,1,1,3,5,4,3,2]
  // const std::string query{"AATGAGAA"};
  // const std::string reference{"AATGA"};
  // const std::string query{"AATGAGAATTTTTTTTTTTTTTTTT"};
  // const std::string reference{"AATGAAAAAAAAAAAAAAAAAA"};
  // int score = xdrop(query, reference, true);

  // int i = 12;
  // const std::string query = TEST_queries[i];
  // const std::string reference = TEST_refs[i];
  // int score = xdrop(query, reference, true);
  // std::cout <<  score << std::endl;

  // Tests  
  std::vector<std::string> refs = TEST_refs;
  std::vector<std::string> qers = TEST_queries;

  for (auto& [path_a, path_b] : INPUT_BATCHS) {
    auto refsnew = loadSequences(path_a);
    auto queriesnew = loadSequences(path_b);
    refs.insert(refs.end(), refsnew.begin(), refsnew.end());
    qers.insert(qers.end(), queriesnew.begin(), queriesnew.end());
  }
  // std::vector<std::string> sequences(refs.size() + qers.size());
  // std::vector<ipu::Comparison> comparisons(refs.size());
  std::vector<std::string> sequences(2*refs.size());
  std::vector<ipu::Comparison> comparisons(1*refs.size());
  for (int i = 0; i < refs.size(); ++i) {
    sequences[2*i] = refs[i];
    sequences[2*i+1] = qers[i];
    comparisons[i] = {
      2*i, 2*i+1
    };
    // break;
  }

  PLOGE << "NUMBER OF COMPARISONS: " << refs.size();

  // std::vector<std::string> refs{};
  // std::vector<std::string> qers{};

  // for (int j =0 ; j<10;j++) {
  //   refs.emplace_back(TEST_refs[10]);
  //   qers.emplace_back(TEST_queries[10]);
  // }




  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -1,
    .gapExtend = -1,
    .matchValue = 1,
    .mismatchValue = -1,
    .ambiguityValue = -1,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {1472 /*Tiles used*/, 1000 /*maxAB*/, 200 , 200 * 20, ipu::VertexType::multixdrop, ipu::Algorithm::fillFirst});

  std::vector<ipu::batchaffine::Batch> batches = driver.create_batches(sequences, comparisons);
  std::vector<ipu::batchaffine::BlockAlignmentResults> results;
  for (auto& batch : batches) {
    PLOGI << "NEW BATCH =======================================================";
    ipu::batchaffine::Job* j = driver.async_submit(&batch);
    assert(batch.cellCount > 0);
    assert(batch.dataCount > 0);
    driver.blocking_join(*j);
    results.push_back(batch.get_result());
    delete j;
  }

  return 0;
}
