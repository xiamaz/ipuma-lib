#include <poplar/AvailableVTypes.h>

#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/Vertex.hpp>
#include <type_traits>
// #include <print.h>

// #include <print.h>
#include "poplar/TileConstants.hpp"

inline int max(int a, int b) {
  return a > b ? a : b;
}

inline int min(int a, int b) {
  return a > b ? b : a;
}

const int GAP_PENALTY = 1;
// const int X = 20;
const int neginf = -9999;

typedef int sType;

template<int X>
class MultiBandXDrop : public poplar::MultiVertex {
 private:
  poplar::Output<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>> K1;
  poplar::Output<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>> K2;

 public:
  // Fields
  poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>> simMatrix;
  poplar::Input<size_t> maxNPerTile;
  poplar::Input<int> gapInit;
  poplar::Input<int> gapExt;
  poplar::Input<int> bufSize;
  poplar::Input<int> maxAB;
  poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Seqs;
  poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Meta;
  poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> score;

  private:
  // This is an eventual atomic counter
  int currentN; 

  public:

  bool compute(unsigned workerId) {

    int gI = *gapInit;
    int gE = *gapExt;

    uint8_t* cSeqs = (uint8_t*)&(Seqs[0]);
    // This is an eventual atomic counter
    int myN = workerId;
    volatile int *globalN = &currentN;
    *globalN = 6;

    for (; myN < maxNPerTile;) {
      const int n = myN;
      int lastNoGap, prevNoGap;

      auto a_len = Meta[4 * n];
      int j_offset = Meta[4 * n + 1];
      auto b_len = Meta[4 * n + 2];
      int i_offset = Meta[4 * n + 3];
      int M = a_len;
      int N = b_len;

      uint8_t* a = cSeqs + j_offset;
      uint8_t* b = cSeqs + i_offset;

      if (a_len == 0 || b_len == 0) break;

      int ma = maxAB;
      sType* k1 = &K1[workerId * (ma+2+2)];
      sType* k2 = &K2[workerId * (ma+2+2)];


      // Algo begin
      int T = 0, T_prime = 0;
      {
        unsigned L = 0, U = 0;

        memset(k1, 0, (ma + 2 + 2) * sizeof(sType));
        memset(k2, 0, (ma + 2 + 2) * sizeof(sType));
        int* t1 = &k1[2];
        int* t2 = &k2[2];
        int L1inc = 0, L2inc = 0;

        auto cell_update = [&](int i, int j, int* _t1, int* __t1, int* _t2, int z, int &lastval, int L) {
          int new_lastval = _t1[z];
          int tscore = max(_t2[z] - GAP_PENALTY, _t2[z- 1] - GAP_PENALTY);
          tscore = max(tscore, lastval + simMatrix[a[i]][b[j]]);

          lastval = new_lastval;
          if (tscore < T - X) {
            tscore = neginf;
          }
          __t1[z] = tscore;
          return tscore;
        };

        int k = 0;
        do {
          // if ((int)U-(int)L >= (int)ma) {
          //   printf("L=%d,U=%d ::span=%d,maxAB=%d, SeqV=%d SeqH=%d\n", L, U, U-L, ma, a_len, b_len);
          //   break;
          // }
          // assert(U-L < klen);

          k = k + 1;
          // int lastval = k1[L-1];
          int *_t2 = t2+(-L+L2inc);
          int *_t1 = t1+(-L+L2inc+L1inc);
          int *__t1 = t1+(-L);

          int lastval = _t1[L-1];
          for (size_t i = L; i < U + 1; i++) {
            auto j = k - i - 1;
            int score = cell_update(i, j, _t1, __t1, _t2, i, lastval, L);
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
      }
      score[n] = T;
      // This operation can cause raise conditions, but that will not produce
      // incurrect results. However, double work can be done, which we accept.
      myN = *globalN;
      (*globalN) += 1;
    }
    return true;
  }
};

template class MultiBandXDrop<10>;
template class MultiBandXDrop<20>;
template class MultiBandXDrop<50>;
template class MultiBandXDrop<100>;