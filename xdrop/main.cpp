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

const std::string basepath = "/global/D1/projects/ipumer/inputs_ab/";

const std::string SEQ_H = basepath + "seqs_seed_H.txt";
const std::string SEED_H = basepath + "seeds_H.txt";
const std::string SEQ_V = basepath + "seqs_seed_V.txt";
const std::string SEED_V = basepath + "seeds_V.txt";

using json = nlohmann::json;

std::vector<size_t> loadSeedpos(const std::string& path) {
  std::ifstream file(path);
  int i;
  std::vector<size_t> seeds;
  while (file >> i) {
    seeds.push_back(i);
  }
  return seeds;
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

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::verbose, &consoleAppender);

  std::vector<swatlib::TickTock> loadTimers(3);
  loadTimers[0].tick();
  PLOGE << "PREPARING COMPARISONS";
  auto [seqs, cmps] = ipu::prepareComparisons(SEQ_H, SEQ_V, SEED_H, SEED_V);
  loadTimers[0].tock();
  PLOGE << "TOOK " << loadTimers[0].duration() / 1000 << " seconds";

  ipu::getDatasetStats(seqs, cmps);

  auto driver = ipu::batchaffine::SWAlgorithm({
    .gapInit = -1,
    .gapExtend = -1,
    .matchValue = 1,
    .mismatchValue = -1,
    .ambiguityValue = -1,
    .similarity = swatlib::Similarity::nucleicAcid,
    .datatype = swatlib::DataType::nucleicAcid,
  }, {
    .numVertices = 1472,
    .maxSequenceLength = 10000,
    .maxComparisonsPerVertex = 100,
    .vertexBufferSize = 10000 * 2 * 6,
    .vtype = ipu::VertexType::multixdrop,
    .fillAlgo = ipu::Algorithm::greedy
  });

  PLOGE << "CREATE BATCHES";
  loadTimers[1].tick();
  auto batches = driver.create_batches(seqs, cmps);
  loadTimers[1].tock();
  PLOGE << "TOOK " << loadTimers[1].duration() / 1000 << " seconds";
  std::vector<ipu::BlockAlignmentResults> results;
  for (auto& batch : batches) {
    ipu::batchaffine::Job* j = driver.async_submit(&batch);
    assert(batch.cellCount > 0);
    assert(batch.dataCount > 0);
    driver.blocking_join(*j);
    results.push_back(batch.get_result());
    delete j;
  }

  return 0;
}
