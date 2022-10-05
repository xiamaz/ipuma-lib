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
 
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

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


seqan3::dna4_vector convertSequence(const std::string& string) {
    seqan3::dna4_vector dna4_str{};
    for (auto c : string) dna4_str.push_back(seqan3::assign_char_to(c, seqan3::dna4{}));
    return dna4_str;
}

std::vector<seqan3::dna4_vector> convertSequences(const std::vector<std::string>& strings) {
    std::vector<seqan3::dna4_vector> seqs{};
    for (const auto& s : strings) {
        seqs.push_back(convertSequence(s));
    }
    return seqs;
}

std::vector<int> seqanAlign(const std::vector<std::string>& queryStrs, const std::vector<std::string>& referenceStrs) {
    auto queries = convertSequences(queryStrs);
    auto references = convertSequences(referenceStrs);
    // seqan3::dna4_vector query = convertSequence(queryStr);
    // seqan3::dna4_vector reference = convertSequence(referenceStr);

    // Configure the alignment kernel.
    auto config =
        seqan3::align_cfg::method_global(
            seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
            seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}
        ) | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-1}};;

    auto results = seqan3::align_pairwise(seqan3::views::zip(queries, references), config);
    std::vector<int> scores{};
    for (auto const& res : results) {
        scores.push_back(res.score());
    }
    return scores;
}


int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

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
  }, {1, 200, 40, 200 * 20, ipu::VertexType::xdrop, ipu::Algorithm::fillFirst});
  driver.compare_local(refs, qers);
  auto aln_results = driver.get_result();

  // driver.compare_local(queries, refs);
  // auto aln_results = driver.get_result();
  // checkResults(aln_results);

  std::cout <<  "We\t" << "We2k\t" << "IPU\t" << "Equal" << endl;
  auto scores_seqan = seqanAlign(refs, qers);
  for (size_t i = 0; i < qers.size(); i++) {
    const std::string query = qers[i];
    const std::string reference = refs[i];
    int score = xdrop(query, reference, true);
    int score2 = xdrop2k(query, reference, true);
    std::cout << score  << "\t"
    << score2 << "\t"
     << aln_results.scores[i] << "\t"
      <<  (score == aln_results.scores[i]) << endl;
  }

  // for (auto& [path_a, path_b] : INPUT_BATCHS) {
  //   auto refs = loadSequences(path_a);
  //   auto queries = loadSequences(path_b);
    
  //   for (size_t i = 0; i < queries.size(); i++) {
  //     const std::string query = queries[i];
  //     const std::string reference = refs[i];
  //     auto score_seqan = seqanAlign(std::vector<std::string>{reference}, std::vector<std::string>{query});
  //     int score = xdrop(query, reference, true);
  //     std::cout << score << "\t" << scores_seqan[i] << "\t" << (score == scores_seqan[i]) << endl;
  //     if (score != scores_seqan[i]) {
  //       PLOGE.printf("ref=%s, quer=%s", reference.c_str(), query.c_str());
  //     }
  //   }
  // }
  return 0;
}
