#include <poplar/AvailableVTypes.h>

#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/Vertex.hpp>
#include <type_traits>

#include "poplar/TileConstants.hpp"

inline int max(int a, int b) {
  return a > b ? a : b;
}

inline int min(int a, int b) {
  return a > b ? b : a;
}

const int GAP_PENALTY = 1;
const int X = 4;
const int neginf = -9999;

class GreedyXDrop : public poplar::Vertex {
 private:
  poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> VTR0;
  poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> VTR1;
  poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> VTT;

 public:
  // Fields
  poplar::Input<size_t> maxNPerTile;
  poplar::Input<int> gapInit;
  poplar::Input<int> gapExt;
  poplar::Input<int> bufSize;
  poplar::Input<int> maxSequenceLength;
  poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Seqs;
  poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Meta;
  poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> score;

  bool compute() {
    const bool cut = true;

    int gI = *gapInit;
    int gE = *gapExt;

    uint8_t* cSeqs = (uint8_t*)&(Seqs[0]);

    for (int n = 0; n < maxNPerTile; ++n) {
      int lastNoGap, prevNoGap;

      auto a_len = Meta[4 * n];
      int j_offset = Meta[4 * n + 1];
      auto b_len = Meta[4 * n + 2];
      int i_offset = Meta[4 * n + 3];

      uint8_t* a = cSeqs + j_offset;
      uint8_t* b = cSeqs + i_offset;

      if (a_len == 0 || b_len == 0) break;


      // Algo begin
      uint8_t *s1 = a;
      uint8_t *s2 = b;
      int N = a_len;
      int M = b_len;


      int mis = -2;
      int mat = 2;
      int X = 10;
      int xdrop_offset = ((X + mat / 2) / (mat - mis)) + 1;

      int inf = 20000;

      auto calcSPrime = [&](int ij, int dis){
          // return floor(
          //     ij * ((double) mat /2) - dis * (mat - mis)
          // );
          return ij * ( mat /2) - dis * (mat - mis);
      };

      int i = 0;
      for (; i < min(M, N) && s1[i] == s2[i]; ++i){ }
      int R_offset = N;

      // int* R0 = (int*) malloc(2 * N * sizeof(int));
      int* R0 = &VTR0[0];
      // int* R1 = (int*) malloc(2 * N * sizeof(int));
      int* R1 = &VTR1[0];
      // int* T = (int*) malloc((M+N+xdrop_offset+1)*sizeof(int));
      int* T = &VTT[0];
      memset(R0, -inf, 2*N*sizeof(int));
      memset(R1, -inf, 2*N*sizeof(int));
      memset(T, 0, (M+N+xdrop_offset+1)*sizeof(int));

      R0[R_offset] = i;
      T[xdrop_offset] = calcSPrime(i + i, 0);

      int Tt = T[xdrop_offset];
      int d = 0;
      int Lmin = 0;
      int Lmax = -inf;
      int L = 0;
      int Umin = inf;
      int Umax = 0;
      int U = 0;

      if (i == N) Lmax = 0;
      if (i == M) Umin = 0;

      for (int d = 1; d < M + N; ++d) {
          int dd = d - xdrop_offset;
          for (int k = L-1; k <= U+1; ++k) {
              int i = -1;
              if (L < k) {
                  i = max(i, R0[k-1+R_offset]+1);
              }
              if (L <= k && k <= U) {
                  i = max(i, R0[k+R_offset]+1);
              }
              if (k < U) {
                  i = max(i, R0[k+1+R_offset]);
              }
              int j = i - k;
              int cs = calcSPrime(i+j, d);
              // if (i >= 0 && j>=0 && i < M && j < N) H(i, j) = H(i, j)+1; // DEBUG
              // std::cout << i << " " << j << " " << cs << std::endl;
              if (i >= 0 && j >= 0 && cs >= T[dd+xdrop_offset] - X) {
                  while (i < M && j < N && s1[i] == s2[j]) {
                      // H(i, j) = H(i, j)+1; // DEBUG
                      ++i, ++j;
                  }
                  R1[k + R_offset] = i;
                  // std::cout << i << " ,j=" << j << " ,M=" << M << " ,N=" << N << " ,cs=" << cs << std::endl;
                  cs = calcSPrime(i+j, d);
                  Tt = max(Tt, cs);
              } else {
                  R1[k + R_offset] = -inf;
              }
          }
          // std::cout << Tt << std::endl;
          T[d+xdrop_offset] = Tt;

          for (int rk = 0; rk < 2*N; ++rk) {
              int i = R1[rk];
              int k = rk - R_offset;
              if (i == N+k) Lmax = max(Lmax, k);
              if (i == M) Lmax =   min(Umin, k);
              if (i > -inf) {
                  Umax = max(Umax, k);
                  Lmin = min(Lmin, k);
              }
          }
          L = max(Lmin, Lmax+2);
          U = min(Umax, Umin-2);

          int* h = R1;
          R1 = R0;
          R0 = h;

          if (L > U + 2) break;
      }

      // free(R0);
      // free(R1);
      // free(T);


      // Algo End

      score[n] = Tt;
    }
    return true;
  }
};