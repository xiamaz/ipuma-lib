#include <poplar/AvailableVTypes.h>

#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/Vertex.hpp>
#include <type_traits>
#include <print.h>

#include "poplar/TileConstants.hpp"

#include "../src/core/xdrop.hpp"

inline int max(int a, int b) {
  return a > b ? a : b;
}

inline int min(int a, int b) {
  return a > b ? b : a;
}

const int GAP_PENALTY = 1;
// const int neginf = -9999;

typedef int sType;
#define GAP_PENALTY 1


template<int X>
class SeedExtendRestrictedXDrop : public poplar::MultiVertex {
 private:
  poplar::Output<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>> K1;
  poplar::Output<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>> K2;

 public:
  // Fields
  poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>> simMatrix;
  poplar::Input<size_t> maxNPerTile;
  poplar::Input<int> seedLength;
  poplar::Input<int> gapInit;
  poplar::Input<int> gapExt;
  poplar::Input<int> bufSize;
  poplar::Input<int> restrictedSize;
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

    for (; myN < maxNPerTile*2;) {
      const int n = myN>>1;
      const int isLeft = myN % 2 == 0;

      auto a_len =   Meta[6 * n];
      int j_offset = Meta[6 * n + 1];
      auto b_len =   Meta[6 * n + 2];
      int i_offset = Meta[6 * n + 3];
      int a_seed_begin = Meta[6 * n + 4];
      int b_seed_begin  = Meta[6 * n + 5];

      int M = a_len;
      int N = b_len;

      uint8_t* a = cSeqs + j_offset;
      uint8_t* b = cSeqs + i_offset;

      if (a_len == 0 || b_len == 0) break;

      // printf("a_seed_begin %d, b_seed_begin %d\n", a_seed_begin, b_seed_begin);

      int klen = restrictedSize+2+2;
      sType* k1 = &K1[workerId * (klen)];
      sType* k2 = &K2[workerId * (klen)];

      int* _k1 = &k1[2];
      int* _k2 = &k2[2];

      sType partscore = 0;
      if (isLeft) {
        memset(k1, 0, (klen) * sizeof(sType));
        memset(k2, 0, (klen) * sizeof(sType));
        partscore = ipumacore::xdrop::xdrop_smart_restricted_extend_left<X, GAP_PENALTY, poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>>&, sType>(
          a, a_len, a_seed_begin,
          b, b_len, b_seed_begin,
          seedLength,
          simMatrix, _k1, _k2, klen
        ) + seedLength;
        // int x = seedLength;
        // printf("partscore left %d %d = %d\n", partscore - x, x, partscore);
      } else {
        // Right hand side
        memset(k1, 0, (klen) * sizeof(sType));
        memset(k2, 0, (klen) * sizeof(sType));
        partscore = ipumacore::xdrop::xdrop_smart_restricted_extend_right<X, GAP_PENALTY, poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>>&, sType>(
          a, a_len, a_seed_begin,
          b, b_len, b_seed_begin,
          seedLength,
          simMatrix, _k1, _k2, klen
        );
        // printf("partscore right %d = %d\n", partscore, partscore);
      }

      score[n] += partscore;
      // This operation can cause raise conditions, but that will not produce
      // incurrect results. However, double work can be done, which we accept.
      myN = *globalN;
      (*globalN) += 1;
    }

    // if (workerId == 0)
    //   for (int i = 0; i < 5; i++) {
    //     printf("Score[%d] = %d\n", i, score[i]);
    //   }
    return true;
  }
};

template class SeedExtendRestrictedXDrop<10>;
// template class SeedExtendRestrictedXDrop<20>;
// template class SeedExtendRestrictedXDrop<50>;
// template class SeedExtendRestrictedXDrop<100>;