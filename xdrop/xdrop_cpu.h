#pragma once

#include <plog/Log.h>

#include <algorithm>
#include <iostream>
#include <span>
#include <string>
#include <vector>
#include <limits>
#include "../src/swatlib/swatlib.h"

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

const int GAP_PENALTY = 1;
const int X = 10;
const int neginf = -9999;

int xdrop(const std::string& query, const std::string& reference, bool cut) {
  auto encoder = getEncoder(DataType::nucleicAcid);
  auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
  auto ref = encoder.encode(reference);
  auto quer = encoder.encode(query);

  int T_prime = 0, T = 0, L = 0, U = 0;

  int M = reference.length();
  int N = query.length();
  Matrix<int> H(M + 1, N + 1, 0); // DEBUG
  Matrix<int> C(M + 1, N + 1, 0); // DEBUG

  // Can also be malloc with: k1[0] = 0, k2[0:2] = 0
  int* k1 = &((int*) calloc(M + 2, sizeof(int)))[1];
  int* k2 = &((int*) calloc(M + 2, sizeof(int)))[1];
  int* k3 = &((int*) calloc(M + 2, sizeof(int)))[1];

  auto cell_update = [&](int i, int j, int* k1, int* k2, int* k3, int z) {
    auto [index, score] = maxtuple({k2[z] - GAP_PENALTY,
                                    k2[z - 1] - GAP_PENALTY,
                                    k1[z - 1] + sim(ref[i], quer[j])});
    if (score < T - X) {
      score = neginf;
    }
    k3[z] = score;
    H(i+1, j+1) = score; // DEBUG
    return score;
  };

  auto rotate = [&]() {
    int* t;  // 1->3, 2->1, 3->2
    t = k1;
    k1 = k2;
    k2 = k3;
    k3 = t;
  };

  int k = 0;
  do {
    k = k + 1;
    for (size_t i = L; i < U + 1; i++) {
      auto j = k - i - 1;
      int score = cell_update(i, j, k1, k2, k3, i);
      C(i + 1, j + 1) = 1; // DEBUG
      T_prime = max(T_prime, score);
    }

    int minL = 99999;
    for (size_t i = L; i < U + 1; i++) {
      int s = k3[i];
      if (s > neginf) {
        minL = i;
        break;
      }
    }

    int maxU = 0;
    for (size_t i = L; i < U + 1; i++) {
      int s = k3[i];
      if (s > neginf) {
        maxU = i;
      }
    }

    if (cut) {
      L = minL;
      U = maxU + 1;
    } else {
      L = 0;
      U = U + 1;
    }

    L = max(L, k + 1 - N);
    U = min(U, M - 1);
    T = T_prime;
    rotate();
  } while (L <= U + 1);

  // PLOGD << H.toString(); // DEBUG
  // PLOGD << C.toString(); // DEBUG

  free(&k3[-1]);
  free(&k2[-1]);
  free(&k1[-1]);
  return T;
}

int xdrop2k(const std::string& query, const std::string& reference, bool cut) {
  auto encoder = getEncoder(DataType::nucleicAcid);
  auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
  auto ref = encoder.encode(reference);
  auto quer = encoder.encode(query);

  int T_prime = 0, T = 0, L = 0, U = 0;

  int M = reference.length();
  int N = query.length();
  Matrix<int> H(M + 1, N + 1, 0); // DEBUG
  Matrix<int> C(M + 1, N + 1, 0); // DEBUG

  // Can also be malloc with: k1[0] = 0, k2[0:2] = 0
  int* k1 = &((int*) calloc(M + 2, sizeof(int)))[1];
  int* k2 = &((int*) calloc(M + 2, sizeof(int)))[1];

  auto cell_update = [&](int i, int j, int* k1, int* k2, int z, int &lastval) {
    int new_lastval = k1[z] + sim(ref[i], quer[j]);
    auto [index, score] = maxtuple({k2[z] - GAP_PENALTY,
                                    k2[z - 1] - GAP_PENALTY,
                                    lastval});
    lastval = new_lastval;
    if (score < T - X) {
      score = neginf;
    }
    k1[z] = score;
    H(i+1, j+1) = score; // DEBUG
    return score;
  };

  auto rotate = [&]() {
    int* t;  // 1->3, 2->1, 3->2
    t = k1;
    k1 = k2;
    k2 = t;
  };

  int k = 0;
  do {
    k = k + 1;
    int lastval = k1[L-1];
    for (size_t i = L; i < U + 1; i++) {
      auto j = k - i - 1;
      int score = cell_update(i, j, k1, k2, i, lastval);
      C(i + 1, j + 1) = 1; // DEBUG
      T_prime = max(T_prime, score);
    }

    int minL = 99999;
    for (size_t i = L; i < U + 1; i++) {
      int s = k1[i];
      if (s > neginf) {
        minL = i;
        break;
      }
    }

    int maxU = 0;
    for (size_t i = L; i < U + 1; i++) {
      int s = k1[i];
      if (s > neginf) {
        maxU = i;
      }
    }

    if (cut) {
      L = minL;
      U = maxU + 1;
    } else {
      L = 0;
      U = U + 1;
    }

    L = max(L, k + 1 - N);
    U = min(U, M - 1);
    T = T_prime;
    rotate();
  } while (L <= U + 1);

  // PLOGD << H.toString(); // DEBUG
  // PLOGD << C.toString(); // DEBUG

  free(&k2[-1]);
  free(&k1[-1]);
  return T;
}