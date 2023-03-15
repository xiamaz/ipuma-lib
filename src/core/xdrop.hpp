#pragma once

#ifndef __POPC__
#include <plog/Log.h>

#include <iostream>

#include "../src/swatlib/swatlib.h"
using namespace swatlib;
#include <string>
#include <vector>
#endif

#include <type_traits>

#ifdef __IPU__
#include <print.h>

#include <Alignment.hpp>
#include <ipu_memory_intrinsics>
#include <ipu_vector_math>
#endif

#include <algorithm>
#include <limits>

template <typename T>
inline T max(T a, T b) {
#ifdef __IPU__
  if constexpr (std::is_same<T, float>()) {
    return __builtin_ipu_max(a, b);
  } else {
    return a > b ? a : b;
  }
#else
  return a > b ? a : b;
#endif
}

template <typename T>
inline T max(T a, T b, T c) {
  return max<T>(max<T>(a, b), c);
}

template <typename T>
inline T min(T a, T b) {
#ifdef __IPU__
  if constexpr (std::is_same<T, float>()) {
    return __builtin_ipu_min(a, b);
  } else {
    return a > b ? b : a;
  }
#else
  return a > b ? b : a;
#endif
}

template <typename T>
inline T min(T a, T b, T c) {
  return min<T>(min<T>(a, c), c);
}


namespace ipumacore {
namespace xdrop {

#ifndef __POPC__
#include "../src/swatlib/swatlib.h"
namespace debug {

template <int X, int GAP_PENALTY>
int xdrop(const std::string& query, const std::string& reference, bool cut) {
  const int neginf = -999999;
  auto encoder = getEncoder(DataType::nucleicAcid);
  auto sim = selectMatrix(Similarity::nucleicAcid, 1, -1, -1);
  auto ref = encoder.encode(reference);
  auto quer = encoder.encode(query);

  int T_prime = 0, T = 0, L = 0, U = 0;

  int M = reference.length();
  int N = query.length();
  Matrix<int> H(M + 1, N + 1, 0);  // DEBUG
  Matrix<int> C(M + 1, N + 1, 0);  // DEBUG

  // Can also be malloc with: k1[0] = 0, k2[0:2] = 0
  int* k1 = &((int*)calloc(M + 2, sizeof(int)))[1];
  int* k2 = &((int*)calloc(M + 2, sizeof(int)))[1];
  int* k3 = &((int*)calloc(M + 2, sizeof(int)))[1];

  auto cell_update = [&](int i, int j, int* k1, int* k2, int* k3, int z) {
    auto score = max(k2[z] - GAP_PENALTY,
                     k2[z - 1] - GAP_PENALTY,
                     k1[z - 1] + sim(ref[i], quer[j]));

#ifdef PRINT_DEBUG
    PLOGW.printf("Q: %s, R: %s", query.c_str(), reference.c_str());
    PLOGW.printf("i %d, j %d, k1[z - 1] (%d) + sim(ref[i] (%c), quer[j] (%c)) (%d) => %d", i, j, k1[z - 1], reference[i], query[j], sim(ref[i], quer[j]), score);
#endif
    if (score < T - X) {
      score = neginf;
    }
    k3[z] = score;
    H(i + 1, j + 1) = score;  // DEBUG
    // PLOGE <<  H.toString();
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
      C(i + 1, j + 1) = 1;  // DEBUG
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
}  // namespace debug
#endif

template <int X, int GAP_PENALTY, bool reversed, typename simT, typename sType>
inline __attribute__((always_inline)) sType xdrop_doubleband(const uint8_t* query, int N, const uint8_t* reference, int M, simT sim, sType* k1, sType* k2) {
  const sType neginf = -999999;
  // printf("NORMAL: N = %d, M = %d\n", N, M);
  if (N == 0 || M == 0) {
    return 0;
  }
  // for (size_t i = 0; i < tmp; i++) {
  //   if (k1[i] != 0) {
  //     printf("k1[%d] == %d\n", i, k1[i]);
  //   }
  //   if (k2[i] != 0) {
  //     printf("k2[%d] == %d\n", i, k2[i]);
  //   }
  // }
  

  int L = 0, U = 0;
  sType T_prime = 0, T = 0;

#define adrREF(zzz) ((!reversed) ? zzz : (M - 1) - zzz)
#define adrQER(zzz) ((!reversed) ? zzz : (N - 1) - zzz)
#ifdef PRINT_DEBUG
  Matrix<int> H(M + 1, N + 1, 0);  // DEBUG
  Matrix<int> C(M + 1, N + 1, 0);  // DEBUG
#endif

  auto rotate = [&]() {
    sType* t;  // 1->3, 2->1, 3->2
    t = k1;
    k1 = k2;
    k2 = t;
  };

  sType lastval = 0;
  int k = 0;
  do {
    k = k + 1;
    if (k != 1) {
      lastval = k1[L - 1];
    }
    for (size_t i = L; i < U + 1; i++) {
      int32_t j = k - i - 1;
      sType new_lastval = k1[i];
      sType score = max<sType>(k2[i] - (sType) GAP_PENALTY,
                               k2[i - 1] - (sType) GAP_PENALTY,
                               lastval + sim[reference[adrREF(i)]][query[adrQER(j)]]);
// #ifdef PRINT_DEBUG
//       printf("\n\n\n");
//       printf("Pos(%d, %d):\n", i + 1, j + 1);
//       printf(" k2[i] - (sType) GAP_PENALTY                           %d - %d = %d\n", k2[i], GAP_PENALTY, k2[i] - (sType) GAP_PENALTY);
//       printf(" k2[i - 1] - (sType) GAP_PENALTY                       %d - %d = %d\n", k2[i - 1], GAP_PENALTY, k2[i - 1] - (sType) GAP_PENALTY);
//       printf(" lastval + sim[reference[adrREF(i)]][query[adrQER(j)]] %d + %d = %d\n", lastval, sim[reference[adrREF(i)]][query[adrQER(j)]], lastval + sim[reference[adrREF(i)]][query[adrQER(j)]]);
// #endif
      lastval = new_lastval;
      if (score < T - X) {
        score = neginf;
      }
      k1[i] = score;
#ifdef PRINT_DEBUG
      H(i + 1, j + 1) = score;  // DEBUG
      C(i + 1, j + 1) = 1;      // DEBUG

  // printf("k1 = ");
  // for (int i = 0; i < std::max(M, N)+2; i++)  {
  //   printf("%7d ", k1[i - 1]);
  // }
  // printf("\n");
  // printf("k2 = ");
  // for (int i = 0; i < std::max(M, N)+2; i++)  {
  //   printf("%7d ", k2[i - 1]);
  // }
  // printf("\n");

  // PLOGD << H.toString();  // DEBUG
  // PLOGD << C.toString();  // DEBUG

#endif
      T_prime = max<sType>(T_prime, score);
    }

    int minL = 99999;
    for (size_t i = L; i < U + 1; i++) {
      sType s = k1[i];
      if (s > neginf) {
        minL = i;
        break;
      }
    }

    int maxU = 0;
    for (size_t i = L; i < U + 1; i++) {
      sType s = k1[i];
      if (s > neginf) {
        maxU = i;
      }
    }

    L = minL;
    U = maxU + 1;

    L = max<int>(L, k + 1 - N);
    U = min<int>(U, M - 1);
    T = T_prime;
    rotate();
  } while (L <= U + 1);

#ifdef PRINT_DEBUG
  PLOGD << H.toString();  // DEBUG
  PLOGD << C.toString();  // DEBUG
#endif

  return T;
}

#ifndef __POPC__
template <int X, int GAP_PENALTY, bool reversed>
int xdrop_doubleband_cpu(const std::vector<uint8_t>& query, const std::vector<uint8_t>& reference, Matrix<int8_t> sim) {
  int M = reference.size();
  int N = query.size();
  int s = std::max(M, N) + 10;
  // Can also be malloc with: k1[0] = 0, k2[0:2] = 0
  int* _k1 = (int*) calloc(s, sizeof(int));
  int* _k2 = (int*) calloc(s, sizeof(int));
  for (int i = 0; i < s; i++) {
    _k1[i] = -9999;
    _k2[i] = -9999;
  }

  int* k1 = &((int*)_k1)[5];
  int* k2 = &((int*)_k2)[5];
  int score = xdrop_doubleband<X, GAP_PENALTY, reversed, std::vector<std::vector<int8_t>>, int>(query.data(), N, reference.data(), M, sim.toVector(), k1, k2);
  free(_k1);
  free(_k2);
  return score;
}
#endif

template <int X, int GAP_PENALTY, bool reversed, typename simT, typename sType>
inline __attribute__((always_inline)) int xdrop_doubleband_restricted(const uint8_t* query, int N, const uint8_t* reference, int M, simT sim, sType* t1, sType* t2, int klen) {
  const sType neginf = -999999;
  // printf("RESTRICTED: N = %d, M = %d\n", N, M);
  // for (size_t i = 0; i < tmp; i++) {
  //   if (t1[i] != 0) {
  //     printf("t1[%d] == %d\n", i, t1[i]);
  //   }
  //   if (t2[i] != 0) {
  //     printf("t2[%d] == %d\n", i, t2[i]);
  //   }
  // }
  if (M == 0 || N == 0) return 0;
  sType T_prime = 0, T = 0;
  int L = 0, U = 0;
  int L1inc = 0, L2inc = 0;

#ifdef PRINT_DEBUG
  Matrix<int> H(M + 1, N + 1, 0);  // DEBUG
  Matrix<int> C(M + 1, N + 1, 0);  // DEBUG
#endif

#define adrREF(zzz) ((!reversed) ? zzz : (M - 1) - zzz)
#define adrQER(zzz) ((!reversed) ? zzz : (N - 1) - zzz)

  sType lastval = 0;
  int k = 0;
  do {
#ifndef __POPC__
    if (U - L >= klen) {
      printf("L=%d,U=%d ::span=%d\n", L, U, U - L);
      printf("SeqV=%d SeqH=%d\n", N, M);
    }
    assert(U - L < klen);
#endif
    k = k + 1;
    // int lastval = k1[L-1];
    sType* _t2 = t2 + (-L + L2inc);
    sType* _t1 = t1 + (-L + L2inc + L1inc);
    sType* __t1 = t1 + (-L);

    if (k != 1) {
      lastval = _t1[L - 1];
    }
    // sType lastval = _t1[L - 1];
    for (size_t i = L; i < U + 1; i++) {
      auto j = k - i - 1;

      sType new_lastval = _t1[i];
      sType score = max<sType>(
          _t2[i] -  (sType) GAP_PENALTY,
          _t2[i - 1] -  (sType) GAP_PENALTY,
          lastval + sim[reference[adrREF(i)]][query[adrQER(j)]]);
      // printf("score[%d]: %d = max(%d, %d, %d)\n", i, score, _t2[i] - GAP_PENALTY, _t2[i - 1] - GAP_PENALTY, lastval + sim[reference[adrREF(i)]][query[adrQER(j)]]);
      // printf("   lv: %d, sim %d,[v %d, i %d][v %d, i %d]\n", lastval,  sim[reference[adrREF(i)]][query[adrQER(j)]], reference[adrREF(i)], adrREF(i), query[adrQER(j)], adrQER(j));

      lastval = new_lastval;
      if (score < T - X) {
        score = neginf;
      }
      __t1[i] = score;

#ifdef PRINT_DEBUG
      H(i + 1, j + 1) = score;  // DEBUG
      C(i + 1, j + 1) = 1;       // DEBUG
#endif

      T_prime = max<sType>(T_prime, score);
    }
    // printf("T_prime[%d]: %d\n",k,  T_prime);
    // if (k == 5) {break;}

    int minL = 99999;
    for (size_t i = L; i < U + 1; i++) {
      sType s = __t1[i];
      if (s > neginf) {
        minL = i;
        break;
      }
    }

    int maxU = 0;
    for (size_t i = L; i < U + 1; i++) {
      sType s = __t1[i];
      if (s > neginf) {
        maxU = i;
      }
    }

    int oldL = L;
    t1[U - L + 1] = neginf;
    L = max(minL, k + 1 - N);
    U = min(maxU + 1, M - 1);
    T = T_prime;
    L1inc = L - oldL;

    // Rotate
    sType* t;  // 1->3, 2->1, 3->2
    t = t1;
    t1 = t2;
    t2 = t;

    int z;
    z = L1inc;
    L1inc = L2inc;
    L2inc = z;
  } while (L <= U + 1);

#ifdef PRINT_DEBUG
  PLOGD << "EXPECTED: " << H.toString();  // DEBUG
  PLOGD << "EXPECTED: " << C.toString();  // DEBUG
#endif
  return T;
}

#ifndef __POPC__
template <int X, int GAP_PENALTY, int klen, bool reversed>
int xdrop_doubleband_restricted_cpu(const std::vector<uint8_t>& query, const std::vector<uint8_t>& reference, Matrix<int8_t> sim) {
  int M = reference.size();
  int N = query.size();
  int s = klen + 10;
  int* _t1 = (int*) calloc(s, sizeof(int));
  int* _t2 = (int*) calloc(s, sizeof(int));
  int *t1 = &_t1[5];
  int *t2 = &_t2[5];

  for (int i = 0; i < s; i++) {
    _t1[i] = -9999;
    _t2[i] = -9999;
  }

  auto score = xdrop_doubleband_restricted<X, GAP_PENALTY, reversed, std::vector<std::vector<int8_t>>, int>(query.data(), N, reference.data(), M, sim.toVector(), t1, t2, klen);
  free(_t1);
  free(_t2);
  return score;
}
#endif

// Double Band
template <int X, int GAP_PENALTY, typename simT, typename sType>
inline __attribute__((always_inline)) sType xdrop_extend_right(const uint8_t* a, int a_len, int a_seed_begin, const uint8_t* b, int b_len, int b_seed_begin, int seedLength, simT sim, sType* k1, sType* k2) {
  return xdrop_doubleband<X, GAP_PENALTY, false, simT, sType>(
      a + seedLength + a_seed_begin, a_len - seedLength - a_seed_begin,
      b + seedLength + b_seed_begin, b_len - seedLength - b_seed_begin,
      sim, k1, k2);
}

template <int X, int GAP_PENALTY, typename simT, typename sType>
inline __attribute__((always_inline)) sType xdrop_extend_left(const uint8_t* a, int a_len, int a_seed_begin, const uint8_t* b, int b_len, int b_seed_begin, int seedLength, simT sim, sType* k1, sType* k2) {
  return xdrop_doubleband<X, GAP_PENALTY, true, simT, sType>(
      a, a_seed_begin,
      b, b_seed_begin,
      sim, k1, k2);
}

// Restricted
template <int X, int GAP_PENALTY, typename simT, typename sType>
inline __attribute__((always_inline)) sType xdrop_restricted_extend_right(const uint8_t* a, int a_len, int a_seed_begin, const uint8_t* b, int b_len, int b_seed_begin, int seedLength, simT sim, sType* k1, sType* k2, int klen) {
  return xdrop_doubleband_restricted<X, GAP_PENALTY, false, simT, sType>(
      a + seedLength + a_seed_begin, a_len - seedLength - a_seed_begin,
      b + seedLength + b_seed_begin, b_len - seedLength - b_seed_begin,
      sim, k1, k2, klen);
}

template <int X, int GAP_PENALTY, typename simT, typename sType>
inline __attribute__((always_inline)) sType xdrop_restricted_extend_left(const uint8_t* a, int a_len, int a_seed_begin, const uint8_t* b, int b_len, int b_seed_begin, int seedLength, simT sim, sType* k1, sType* k2, int klen) {
  return xdrop_doubleband_restricted<X, GAP_PENALTY, true, simT, sType>(
      a, a_seed_begin,
      b, b_seed_begin,
      sim, k1, k2, klen);
}

// Smart Switching
template <int X, int GAP_PENALTY, typename simT, typename sType>
inline __attribute__((always_inline)) sType xdrop_smart_restricted_extend_left(const uint8_t* a, int a_len, int a_seed_begin, const uint8_t* b, int b_len, int b_seed_begin, int seedLength, simT sim, sType* k1, sType* k2, int klen) {
  const int a_len_real = a_seed_begin;
  const int b_len_real = b_seed_begin;
  if (a_seed_begin < b_seed_begin) {
    std::swap(a, b);
    std::swap(a_seed_begin, b_seed_begin);
  }
  int mlen = min<int>(a_len_real, b_len_real);

  if (mlen >= klen) {
    return xdrop_doubleband_restricted<X, GAP_PENALTY, true, simT, sType>(
        a, a_seed_begin,
        b, b_seed_begin,
        sim, k1, k2, klen);

  } else {
    return xdrop_doubleband<X, GAP_PENALTY, true, simT, sType>(
        a, a_seed_begin,
        b, b_seed_begin,
        sim, k1, k2);
  }
}

template <int X, int GAP_PENALTY, typename simT, typename sType>
inline __attribute__((always_inline)) sType xdrop_smart_restricted_extend_right(const uint8_t* a, int a_len, int a_seed_begin, const uint8_t* b, int b_len, int b_seed_begin, int seedLength, simT sim, sType* k1, sType* k2, int klen) {
  const int a_len_real = a_len - seedLength - a_seed_begin;
  const int b_len_real = b_len - seedLength - b_seed_begin;
  if (a_len_real < b_len_real) {
    std::swap(a, b);
    std::swap(a_seed_begin, b_seed_begin);
    std::swap(a_len, b_len);
  }
  int mlen = min<int>(a_len_real, b_len_real);
  if (mlen >= klen) {
    return xdrop_doubleband_restricted<X, GAP_PENALTY, false, simT, sType>(
        a + seedLength + a_seed_begin, a_len - seedLength - a_seed_begin,
        b + seedLength + b_seed_begin, b_len - seedLength - b_seed_begin,
        sim, k1, k2, klen);
  } else {
    return xdrop_doubleband<X, GAP_PENALTY, false, simT, sType>(
        a + seedLength + a_seed_begin, a_len - seedLength - a_seed_begin,
        b + seedLength + b_seed_begin, b_len - seedLength - b_seed_begin,
        sim, k1, k2);
  }
}

#ifndef __POPC__
namespace debug {
template <int X, int GAP_PENALTY, int seedLength>
int seed_extend_cpu(const std::vector<uint8_t>& query, int querySeedBeginPos, const std::vector<uint8_t>& reference, int referenceSeedBeginPos, Matrix<int8_t> sim) {
  auto decoder = swatlib::getEncoder(DataType::nucleicAcid);
  PLOGW << decoder.decode(query);
  PLOGW << decoder.decode(reference);
  PLOGD << querySeedBeginPos << ", " << referenceSeedBeginPos;
  auto leftH = std::vector<uint8_t>(query.begin(), query.begin() + querySeedBeginPos);
  auto leftV = std::vector<uint8_t>(reference.begin(), reference.begin() + referenceSeedBeginPos);

  PLOGW << decoder.decode(leftH);
  PLOGW << decoder.decode(leftV);
  auto scoreLeft = xdrop_doubleband_cpu<X, GAP_PENALTY, true>(
      leftH,
      leftV,
      sim);

  auto rightH = std::vector<uint8_t>(query.begin() + querySeedBeginPos + seedLength, query.end());
  auto rightV = std::vector<uint8_t>(reference.begin() + referenceSeedBeginPos + seedLength, reference.end());
  PLOGW << decoder.decode(rightH);
  PLOGW << decoder.decode(rightV);
  auto scoreRight = xdrop_doubleband_cpu<X, GAP_PENALTY, false>(
      rightH,
      rightV,
      sim);
  PLOGE << scoreLeft << "," << scoreRight;
  return scoreLeft + scoreRight + seedLength;
}
}  // namespace debug

template <int X, int GAP_PENALTY>
int seed_extend_cpu(const std::vector<uint8_t>& query, int querySeedBeginPos, const std::vector<uint8_t>& reference, int referenceSeedBeginPos, int seedLength, Matrix<int8_t> sim) {
  int M = reference.size();
  int N = query.size();
  assert(M >= N);
  // Can also be malloc with: k1[0] = 0, k2[0:2] = 0
  int* k1 = &((int*)malloc((M + 2) * sizeof(int)))[0];
  int* k2 = &((int*)malloc((M + 2) * sizeof(int)))[0];
  int* _k1 = &k1[1];
  int* _k2 = &k2[1];

  memset(k1, 0, sizeof(int) * M + 2);
  memset(k2, 0, sizeof(int) * M + 2);
  auto score_right = xdrop_extend_right<X, GAP_PENALTY, std::vector<std::vector<int8_t>>, int>(query.data(), N, querySeedBeginPos, reference.data(), M, referenceSeedBeginPos, seedLength, sim.toVector(), _k1, _k2);

  memset(k1, 0, sizeof(int) * M + 2);
  memset(k2, 0, sizeof(int) * M + 2);
  auto score_left = xdrop_extend_left<X, GAP_PENALTY, std::vector<std::vector<int8_t>>, int>(query.data(), N, querySeedBeginPos, reference.data(), M, referenceSeedBeginPos, seedLength, sim.toVector(), _k1, _k2);

  free(k2);
  free(k1);

  return score_right + score_left + seedLength;
}

template <int X, int GAP_PENALTY>
int seed_extend_restricted_cpu(const std::vector<uint8_t>& query, int querySeedBeginPos, const std::vector<uint8_t>& reference, int referenceSeedBeginPos, int seedLength, Matrix<int8_t> sim, int klen) {
  int N = query.size();
  int M = reference.size();
  // printf("%d => (%d|%d)\n",query.size(),  querySeedBeginPos, query.size() - seedLength - querySeedBeginPos);
  // printf("%d => (%d|%d)\n", reference.size(), referenceSeedBeginPos, query.size() - seedLength - referenceSeedBeginPos);

  // assert(M >= N);
  // Can also be malloc with: k1[0] = 0, k2[0:2] = 0
  size_t buflen = klen + 4 + 4;
  int* k1 = &((int*)malloc((buflen) * sizeof(int)))[0];
  int* k2 = &((int*)malloc((buflen) * sizeof(int)))[0];
  int* _k1 = &k1[4];
  int* _k2 = &k2[4];

  memset(k1, 0, sizeof(int) * buflen);
  memset(k2, 0, sizeof(int) * buflen);
  auto score_left = xdrop_smart_restricted_extend_left<X, GAP_PENALTY, std::vector<std::vector<int8_t>>, int>(
      query.data(), N, querySeedBeginPos,
      reference.data(), M, referenceSeedBeginPos,
      seedLength, sim.toVector(), _k1, _k2, klen);

  memset(k1, 0, sizeof(int) * buflen);
  memset(k2, 0, sizeof(int) * buflen);
  auto score_right = xdrop_smart_restricted_extend_right<X, GAP_PENALTY, std::vector<std::vector<int8_t>>, int>(
      query.data(), N, querySeedBeginPos,
      reference.data(), M, referenceSeedBeginPos,
      seedLength, sim.toVector(), _k1, _k2, klen);

  free(&k2[0]);
  free(&k1[0]);

  // std::cout << score_left << " " << seedLength << " " << score_right << std::endl;
  return score_right + score_left + seedLength;
}
#endif

}  // namespace xdrop
}  // namespace ipumacore