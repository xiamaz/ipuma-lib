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
  // const std::string query{"AATGAGAA"};
  // const std::string reference{"AATGA"};
  const std::string query = queries[5];
  const std::string reference = refs[5];




  int X = 4;
  int neginf = -99;
  int T_prime = 0, T = 0, L = 0, U = 0; 

  int m = reference.length() + 1;
  int n = query.length() + 1;
  Matrix<int> H(m, n, 0);
  Matrix<int> C(m, n, 0);

  int * k1 = &((int*) calloc(min(m, n) + 1, sizeof(int)))[1];
  int * k2 = &((int*) calloc(min(m, n) + 1, sizeof(int)))[1];
  int * k3 = &((int*) calloc(min(m, n) + 1, sizeof(int)))[1];

  auto cell_update = [&](int i, int j, int* k1, int* k2, int* k3, int z){
    auto [index, score] = maxtuple({
      k2[z] - GAP_PENALTY,
      k2[z-1] - GAP_PENALTY,
      k1[z-1] + simpleSimilarity(reference[i-1], query[j-1])
    });
    if (score < T - X) {
      // PLOGW.printf("(%d, %d) DROP", i, j);
      score = neginf;
    }
    k3[z] = score;
    H(i, j) = score;
    return score;
  };

  auto rotate = [&]() {
    int *t; // 1->3, 2->1, 3->2
    t = k1;
    k1 = k2;
    k2 = k3;
    k3 = t;
  };

  int M = reference.length();
  int N = query.length();

  int k = 0;
  int c = 0;
  do {
    k = k + 1;
    int lastU = 0;
    for (size_t i = L; i < U+1; i++) {
      auto j = k - i - 1;
      int _j = j + 1;
      int _i = i + 1;
      int score = cell_update(_i, _j, k1, k2, k3, i);
      if (score == neginf) {
        if ((k3[i-1] == neginf && i != 0) || i == 0) {
          L = i;
        }
        if (k3[i-1] != neginf) {
          lastU = i;
        }
      }
      C(_i, _j) = 1;
      T_prime = max(T_prime, score);
    }
    // find first non neg score
    // L = ;
    if (k3[U] == neginf) {
      U = lastU+1;
    }
    PLOGE.printf("L(%d), U(%d), lastU(%d), k3[L](%d), k3[U](%d)", L, U, lastU, k3[L], k3[U]);

    // int minL = max(M, N);
    // for (int i = U; i >= 0; i--) {
    //   int s = H(i, k-i);
    //   if (s > neginf) {
    //     minL = i;
    //   }
    // }

    // int maxU = 0;
    // for (int i = 0; i < U+1; i++) {
    //   int s = H(i, k-i);
    //   if (s > neginf) {
    //     maxU = i;
    //   }
    // }
    PLOGW.printf("L(%d), U(%d)",  L, U);
    L = minL;
    U = maxU;
    L = max(L, k+1-N);
    U = min(U, M-1);
    T = T_prime;
    rotate();
  } while  (L <= U+1);

  PLOGI << H.toString();
  PLOGI << C.toString();

  return 0;
}
