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
const int X = 20;
const int neginf = -9;

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
  auto klen = 90;
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

  int* t1 = &((int*) calloc(klen+5, sizeof(int)))[1];
  int* t2 = &((int*) calloc(klen+5, sizeof(int)))[1];
  int L1inc = 0, L2inc = 0;

  auto cell_update = [&](int i, int j, int* k1, int* k2, int z, int &lastval, int L) {
    // int new_lastval = k1[z];
    int tnew_lastval = t1[(z-L)+L2inc+L1inc];
    // printf("lv=%d[%d], tlv=%d[%d]\n", new_lastval, z,  tnew_lastval, z-L);
    // assert(new_lastval == tnew_lastval);
    // assert(k2[z] == t2[(z-L)+L2inc]);
    // printf("k2[z - 1] = %d, t2[(z-L)+L2inc - 1] = %d\n", k2[z - 1] , t2[(z-L)+L2inc - 1]);
    // assert(k2[z - 1] == t2[(z-L)+L2inc - 1]);
    auto [tindex, tscore] = maxtuple({t2[(z-L)+L2inc] - GAP_PENALTY,
                                    t2[(z-L)+L2inc - 1] - GAP_PENALTY,
                                    lastval + sim(ref[i], quer[j])});

    // auto [index, score] = maxtuple({k2[z] - GAP_PENALTY,
    //                                 k2[z - 1] - GAP_PENALTY,
    //                                 lastval + sim(ref[i], quer[j])});

    // assert(score == tscore);
    lastval = tnew_lastval;
    if (tscore < T - X) {
      tscore = neginf;
    }
    t1[z-L] = tscore;
    k1[z] = tscore;
    H(i+1, j+1) = tscore; // DEBUG
    return tscore;
  };

  auto rotate = [&]() {
    int* t;  // 1->3, 2->1, 3->2
    t = k1;
    k1 = k2;
    k2 = t;

    // For ts
    t = t1;
    t1 = t2;
    t2 = t;
    
    int z;
    z = L1inc;
    L1inc = L2inc;
    L2inc = z;
  };

  int k = 0;
  do {
    if (U-L >= klen) {
      printf("L=%d,U=%d ::span=%d\n", L, U, U-L);
      printf("SeqV=%d SeqH=%d\n", query.length(), reference.length());
    }
    assert(U-L < klen);

    k = k + 1;
    // int lastval = k1[L-1];
    int lastval = t1[((L-1)-L)+L2inc+L1inc];
    for (size_t i = L; i < U + 1; i++) {
      auto j = k - i - 1;
      int score = cell_update(i, j, k1, k2, i, lastval, L);
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

    int tminL = 99999;
    for (size_t i = L; i < U + 1; i++) {
      int s = t1[i-L];
      if (s > neginf) {
        tminL = i;
        break;
      }
    }
    assert(tminL == minL);

    int tmaxU = 0;
    for (size_t i = L; i < U + 1; i++) {
      int s = t1[i-L];
      if (s > neginf) {
        tmaxU = i;
      }
    }

    int maxU = 0;
    for (size_t i = L; i < U + 1; i++) {
      int s = k1[i];
      if (s > neginf) {
        maxU = i;
      }
    }
    assert(maxU == tmaxU);

    int oldL = L;
    int oldU = U;

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
    // printf("k %2d:\t", k);
    // for (size_t i = 0; i < M ; i++) {
    //   if (i == oldL || i == oldU+1) {
    //     printf("|");
    //   } else {
    //     printf(" ");
    //   }
    //   printf("%2d", k1[i]);
    // }

    // printf("\t L=%2d, U=%2d\n", oldL, oldU+1);

    // int off = max(0, oldU-klen);
    // int *_k1 = &k1[-off];

    // printf("View:\t");
    // int *_k1 = &k1[oldL];
    // for (size_t i = 0; i < klen; i++) {
    //   printf(" ");
    //   printf("%2d", _k1[i]);
    // }
    // printf("\n");

    t1[oldU - oldL+1] = 0;
    // printf("T1 :\t");
    // for (size_t i = 0; i < klen; i++) {
    //   printf(" ");
    //   printf("%2d", t1[i]);
    // }
    // printf("\n");

    L1inc = L - oldL;
    // printf("Linc=%2d\n", L1inc);

    // printf("\n");
    // printf("\n");
    
    rotate();
    // printf("T1 :\t");
    // for (size_t i = 0; i < klen; i++) {
    //   printf(" ");
    //   printf("%2d", t1[i]);
    // } printf("\n");

    // printf("T2 :\t");
    // for (size_t i = 0; i < klen; i++) {
    //   printf(" ");
    //   printf("%2d", t2[i]);
    // } printf("\n");


  } while (L <= U + 1);

  // PLOGD << "EXPECTED: " << H.toString(); // DEBUG
  // PLOGD << "EXPECTED: " << C.toString(); // DEBUG

  // free(&k2[-1]);
  // free(&k1[-1]);
  return T;
}

int xdrop3k(const std::string& query, const std::string& reference, bool cut) {
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
    int new_lastval = k1[z];
    auto [index, score] = maxtuple({k2[z] - GAP_PENALTY,
                                    k2[z - 1] - GAP_PENALTY,
                                    lastval + sim(ref[i], quer[j])});
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

  PLOGD << H.toString(); // DEBUG
  PLOGD << C.toString(); // DEBUG

  free(&k2[-1]);
  free(&k1[-1]);
  return T;
}

int xdrop4k(const std::string& query, const std::string& reference, bool cut) {
  auto klen = 90;
  auto encoder = getEncoder(DataType::nucleicAcid);
  auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
  auto ref = encoder.encode(reference);
  auto quer = encoder.encode(query);

  int T_prime = 0, T = 0, L = 0, U = 0;

  int M = reference.length();
  int N = query.length();
  Matrix<int> H(M + 1, N + 1, 0); // DEBUG
  Matrix<int> C(M + 1, N + 1, 0); // DEBUG

  int* t1 = &((int*) calloc(klen+5, sizeof(int)))[1];
  int* t2 = &((int*) calloc(klen+5, sizeof(int)))[1];
  int L1inc = 0, L2inc = 0;

  auto cell_update = [&](int i, int j, int* _t1, int* __t1, int* _t2, int z, int &lastval, int L) {
    int new_lastval = _t1[z];
    int tscore = max(_t2[z] - GAP_PENALTY, _t2[z- 1] - GAP_PENALTY);
    tscore = max(tscore, lastval + sim(ref[i], quer[j]));

    lastval = new_lastval;
    if (tscore < T - X) {
      tscore = neginf;
    }
    __t1[z] = tscore;
    H(i+1, j+1) = tscore; // DEBUG
    return tscore;
  };

  int k = 0;
  #pragma unroll(2)
  do {
    if (U-L >= klen) {
      printf("L=%d,U=%d ::span=%d\n", L, U, U-L);
      printf("SeqV=%d SeqH=%d\n", query.length(), reference.length());
    }
    assert(U-L < klen);

    k = k + 1;
    // int lastval = k1[L-1];
    int *_t2 = t2+(-L+L2inc);
    int *_t1 = t1+(-L+L2inc+L1inc);
    int *__t1 = t1+(-L);

    int lastval = _t1[L-1];
    for (size_t i = L; i < U + 1; i++) {
      auto j = k - i - 1;
      int score = cell_update(i, j, _t1, __t1, _t2, i, lastval, L);
      C(i + 1, j + 1) = 1; // DEBUG
      T_prime = max(T_prime, score);
    }

    int minL = 99999;
    for (size_t i = L; i < U + 1; i++) {
      int s = __t1[i];
      if (s > neginf) {
        minL = i;
        break;
      }
    }

    int maxU = 0;
    for (size_t i = L; i < U + 1; i++) {
      int s = __t1[i];
      if (s > neginf) {
        maxU = i;
      }
    }

    int oldL = L;
    t1[U - L+1] = 0;
    L = max(minL, k + 1 - N);
    U = min(maxU + 1, M - 1);
    T = T_prime;
    L1inc = L - oldL;

    // Rotate
    int* t;  // 1->3, 2->1, 3->2
    t = t1;
    t1 = t2;
    t2 = t;
    
    int z;
    z = L1inc;
    L1inc = L2inc;
    L2inc = z;
  } while (L <= U + 1);

  PLOGD << "EXPECTED: " << H.toString(); // DEBUG
  PLOGD << "EXPECTED: " << C.toString(); // DEBUG

  // free(&k2[-1]);
  // free(&k1[-1]);
  return T;
}