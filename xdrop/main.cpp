#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>

#include <cxxopts.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <span>

#include "swatlib/swatlib.h"

using json = nlohmann::json;
std::vector<std::string> queries = {
    "AATGAGAATGATGTCGTTCGAAATTTGACCAGTCAAACCGCGGGCAATAAGGTCTTCGTTCAGGGCATAGACCTTAATGGGGGCATTACGCAGACTTTCA",
    "ATCTGGCAGGTAAAGATGAGCTCAACAAAGTGATCCAGCATTTTGGCAAAGGAGGCTTTGATGTGATTACTCGCGGTCAGGTGCCACCTAACCCGTCTGA",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGGTATTGGGAGCACG",
    "CCCCGCACCCGCAAGCCGCCGAGAAAAAAAGGATGAGGGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGAACA",
    "CGTGCTCCCAATACCCGCCAGGAATTAGCGCGCTTACGCTCAGAGCTGGCGGAAAAAGAGAGCCAGCAACTGGATGCGTATTTGCAGGCCTTGCGTAATC",
    "CCCCGCACCCGCAAGCCGCCGAGAAAAAAAGGATGAGGGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGAACA",
    "ATATCATCACTCCGATGGACGTTTCGATTATCGGCTGCGTGGGGAATGGCCCAGGTGAGGCGCTGGTTTCTACACTCGGCGTCACCGGCGGGACGAAACA",
    "CCGGCGGTGGGCGCGTCCGCCAGTGCCGGCGCGAGCAGGACGGCGTAGAGCCGGATGACGTGATCCTGCCGCCGGAGAAGGCAAAAGCGAAAGGGCTGAC",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGGTATTGGGAGCACG",
    "TCCGGCTGGCAGAACTTGACCAGTGCCGATAAAGAAAAGATGTGGCAGGCGCAGGGGCGAATGCTCACCGCACAAAGCCTGAAGATTAATGCGTTGCTGC",
    "CCGCCCCCGCGCACACGGTGCGGCCTGTCCCGCGTATACTCGACCAGCGGCGTCCCGCCCAGCTTCATTCCCGCCAGGTAACCGCTGCCATACGTCAGCC",
    "AAACCATTGCGAGGAAGTGGTTCTACTTGCTGTCGCCGCGGGAGAACAGGTTGTGCCTGCCTCTGAACTTGCTGCCGCCATGAAGCAGATTAAAGAACTC",
    "CCTTCCCCCTAACTTTCCGCCCGCCATGAAGCAGATAAAAGAACTCCAGCGCCTGCTCGGAAAGAAAACGATGGAAAATGAACTCCTCAAAGAAGCCGTT",
    "AGATGTGCCGGTCATTAAGCATAAAGCCGATGGTTTCTCCCCGCACTTGCCGCCAGTGACGCCACGGCCAGTCAGAGAAGATCATAACAACCGCTCCAGT",
    "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCCGTAGTTCTTCACCCGTTTTGAATTCGGGTTTGTAT",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCACTTGATGGCAGATAACGCGCC",
};
std::vector<std::string> refs = {
    "AATGAGAATGATGTCNTTCNAAATTTGACCAGTCAAACCGCGGGCAATAAGGTCTTCGTTCAGGGCATAGACCTTAATGGGGGCATTACGCAGACTTTCA",
    "ATCTGGCAGGTAAAGATGAGCTCAACAAAGTGATCCAGCATTTTGGCAAAGGAGGCTTTGATGTGATTACTCGCGGTCAGGTGCCACCTAANNCGTCTGA",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGTTATTGGCAGACAG",
    "GCACCGTCCAGCCAACCGCCGAGAAGAAAAGAATGAGTGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGACAC",
    "CTGTCTGCCAATAACCGCCAGGAATTAGCGCGCTTACGCTCAGAGCTGGCGGAAAAAGAGAGCCAGCAACTGGATGCGTATTTGCAGGCCTTGCGTAATC",
    "GCACCGTCCAGCCAACCGCCGAGAAGAAAAGAATGAGTGCGATACGGATCAGGATATCTACGGTTTTCTGCCCCGCGCCGTTTTGCAGCCAGTTCCAGAA",
    "AATAATAATAATGTCGCAGTCGTCTTCCATGTCATGCCCCAGATATCCAGAACACAACACCCTAACATAGCGTTACTTAAGGGAAATTGACCGCCGACAC",
    "TGGTTTCTACACTCGGCGTCACCGGCGGCAACAAGAA",
    "AGCGCCGGGCGCGCTTCCGCCAGTGCCTGCGCGAGCAGGACGGCGTAGAGCCGGATGACGTGATCCTGCCGCCGGAGAAGGCAAAAGCGAAAGGGCTGAC",
    "GATTACGCAAGGCCTGCAAATACGCATCCAGTTGCTGGCTCTCTTTTTCCGCCAGCTCTGAGCGTAAGCGCGCTAATTCCTGGCGGTTATTGGCAGACAG",
    "TTCGCCGCGCAGAACCTGACCAGTGCCGATAACGAAAAGATGTGGCAGGCGCAGGGGCGAATGCTCACCGCACAAAGCCTGAAGATTAATGCGTTGCTGC",
    "CTGCGCACCGTCTCACGGTGCAGCCTGTCCCGCGTATACTCGACCAGCGGCGTCCCGCCCAGCTTCATTCCCGCCAGGTAACCGCTGCCATACGTCAGC",
    "TAAGCAATACCAGGAAGGAAGTCTTACTGCTGTCGCCGCCGGAGAACAGGTTGTTCCTGCCTCTGAACTTGCTGCCGCCATGAAGCAGATTAAAGAACTC",
    "TCCTGCCTCTGAACTTGCTGCCGCCATGAAGCAGATTAAAGAACTCCAGCGCCTGCTCGGCAAGAAAACGATGGAAAATGAACTCCTCAAAGAAGCCGTT",
    "TGAGTTGCTCGTCATTAAGACGTAAGGCGATGGTTTCTCCCCGCACTTGCCGCCAGTG",
    "CATCGCCCGATTTTCACGTTCGAGAGCGGCGGAGCGGATCGCTCCTTGTTCTTTTTGCCAGGCCAGTAGTTCTTCACCCGTTTTGAATGCGGGTTTGATA",
    "GCCAGGCAAAATCGGCGTTTCTGGCGGCGATGAGCCATGAGATCCGCACACCGCTGTACGGTATTCTCGGCACTGCTCAACTGCTGGCAGATAACCCCGC",
};

using namespace std;
using namespace swatlib;

template <typename T>
inline std::tuple<int, T> maxtuple(std::initializer_list<T> l) {
  T value;
  int index = -1;
  for (auto it = l.begin(); it != l.end(); ++it) {
    if (index == -1 || *it > value) {
      value = *it;
      index = std::distance(l.begin(), it);
    }
  }
  return {index, value};
}

template <typename T>
inline T maxval(std::initializer_list<T> l) {
  auto [i, v] = maxtuple(l);
  return v;
}

template<typename T>
inline T vmax(std::vector<T> vec) {
  return *std::max_element(vec.begin(), vec.end());
}

template<typename T>
inline T vmin(std::vector<T> vec) {
  return *std::min_element(vec.begin(), vec.end());
}

const int GAP_PENALTY = 1;

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);
  auto similarityMatrix = swatlib::selectMatrix(swatlib::Similarity::nucleicAcid, 1, -1, -1);

  //    0 1 2 3 4 5 6 7 n
  // 0 [0,0,0,0,0,0,0,0,0]
  // 1 [0,1,1,0,0,1,0,1,1]
  // 2 [0,1,2,1,0,1,0,1,2]
  // 3 [0,0,1,3,2,1,0,0,1]
  // 4 [0,0,0,2,4,3,2,1,0]
  // m [0,1,1,1,3,5,4,3,2]
  const std::string query{"AATGAGAA"};
  const std::string reference{"AATGA"};
  int X = 3;
  int neginf = -9;
  int T_prime = 0, T = 0, L = 0, U = reference.length();

  // const auto I = 8;
  // const std::string query = queries[I];
  // const std::string reference = refs[I];

  int m = reference.length() + 1;
  int n = query.length() + 1;
  Matrix<int> H(m, n);
  Matrix<int> C(m, n);

  auto k1 = std::vector<int>(min(m, n), 0);
  auto k2 = std::vector<int>(min(m, n), 0);
  auto k3 = std::vector<int>(min(m, n), 0);

  // std::span<int> k4{k1.data() + 1, 2};
  // k4[0] = 99;

  // PLOGD << printVectorD(k1);
  // return 0;

  auto XDropUpdate = [&](){
    T_prime = max(T_prime, vmax(k3));
    for (size_t i = L; i < U; i++) {
      if (k3[i] > neginf) {
        L = i;
        break;
      }
    }
    for (size_t i = U; i < L; i--) {
      if (k3[i] > neginf) {
        U = i;
        break;
      }
    }
    T = T_prime;
    PLOGD.printf("[L, U](%d, %d), T = %d", L, U, T);
  };

  auto cell_update_top = [&](int i, int j, int z) {
    auto [index, score] = maxtuple({
        0,
        k2[z] - GAP_PENALTY,                                     // Left
        k2[z - 1] - GAP_PENALTY,                                 // Top
        k1[z - 1] + (reference[i - 1] == query[j - 1] ? 1 : -1)  // Diag
    });
    if (score < T - X) {
      score = neginf;
      C(i,j) = 1;
    }
    k3[z] = score;
    H(i,j) = score;
  };

  auto cell_update_bottom = [&](int i, int j, int z) {
    auto [index, score] = maxtuple({
        0,
        k2[z + 1] - GAP_PENALTY,                                 // Left
        k2[z] - GAP_PENALTY,                                     // Top
        k1[z + 1] + (reference[i - 1] == query[j - 1] ? 1 : -1)  // Diag
    });
    if (score < T - X) {
      score = neginf;
      C(i,j) = 1;
    }
    k3[z] = score;
    H(i,j) = score;
  };

  auto rotate = [&]() {
    // 1->3, 2->1, 3->2
    k1.swap(k3);
    k2.swap(k1);
  };

  // Split computatioin into: (◸⟋◿) -> (upper antidiag triangle) + (ext. antidiag) + (lower antidiag triangle)
  // Upper antidiagonal matrix
  for (int i = 0; i < m; i++) {
    for (int j = 1; j < i; j++) {
      cell_update_top(j, i - j, j);
    }
    XDropUpdate();
    PLOGI << printVectorD(k3);
    rotate();
  }

  // Band
  for (int i = 0; i < n - m; i++) {
    for (int j = 1; j < m; j++) {
      cell_update_top(j, m + i - j, j);
    }
    XDropUpdate();
    PLOGI << printVectorD(k3);
    rotate();
  }

  // We are still top centered, but k1 is not... this is an edge case
  // TODO: make bottom pinned can fix this
  k1.insert(k1.begin(), 0);
  // Lower antidiagonal matrix
  for (int i = 0; i < (m - 1); i++) {
    for (int j = 1; j < m - i; j++) {
      cell_update_bottom(j + i, n - j, j - 1);
    }
    XDropUpdate();
    PLOGI << printVectorD(k3);
    rotate();
  }

  PLOGI << H.toString();
  PLOGI << C.toString();

  return 0;
}
