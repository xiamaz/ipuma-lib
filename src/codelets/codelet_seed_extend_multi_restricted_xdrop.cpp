#include <ipu_builtins.h>
#include <poplar/AvailableVTypes.h>
#include <print.h>

#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/Vertex.hpp>
#include <type_traits>

#include "../src/core/xdrop.hpp"
#include "../src/shared_types.h"
#include "poplar/TileConstants.hpp"

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

class TileTelemetry {
 private:
  unsigned cycle_counter_l;
  unsigned cycle_counter_u;

 public:
  TileTelemetry() {
    reset();
  }

  unsigned cycle_counter_diff() {
    unsigned now_l = __builtin_ipu_get_scount_l();
    unsigned now_u = __builtin_ipu_get_scount_u();

    unsigned diff_u = now_u - cycle_counter_u;
    unsigned diff_l = now_l - cycle_counter_l;

    if (now_u < cycle_counter_u)
      return __UINT32_MAX__;
    else if (diff_u == 0)
      return diff_l;
    else if (diff_u == 1) {
      if (now_l >= cycle_counter_l)
        return __UINT32_MAX__;
      else {
        return __UINT32_MAX__ - (cycle_counter_l - now_l);
      }
    } else if (diff_u > 1) {
      return __UINT32_MAX__;
    }
    return 0;
  }

  void reset() {
    cycle_counter_l = __builtin_ipu_get_scount_l();
    cycle_counter_u = __builtin_ipu_get_scount_u();
  }
};

constexpr unsigned numberOfBits(unsigned x) {
  return x < 2 ? x : 1 + numberOfBits(x >> 1);
}

template <int X>
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
  poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> ARange;
  poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> BRange;

 private:
  // This is an eventual atomic counter
  int currentN;

 public:
  bool compute(unsigned workerId) {
    // auto tt = TileTelemetry();
    int gI = *gapInit;
    int gE = *gapExt;

    uint8_t* cSeqs = (uint8_t*)&(Seqs[0]);
    // This is an eventual atomic counter
    int myN = workerId;
    volatile int* globalN = &currentN;
    *globalN = 6;

    if (workerId == 0) {
      memset((char*)&score[0], 0, maxNPerTile * sizeof(score[0]));
    }
    const ipu::XDropMeta* metas = reinterpret_cast<const ipu::XDropMeta*>(&Meta[0]);
    for (; myN < maxNPerTile * 2;) {
      const bool isLeft = myN % 2 == 0;
      const int n = myN >> 1;

      const ipu::XDropMeta meta = metas[n];
      const auto a_len = meta.sizeA;
      const auto j_offset = meta.offsetA;
      const auto b_len = meta.sizeB;
      const auto i_offset = meta.offsetA;

      const int32_t a_seed_begin = meta.seedAStartPos;
      const int32_t b_seed_begin = meta.seedBStartPos;

      const uint8_t* a = cSeqs + j_offset;
      const uint8_t* b = cSeqs + i_offset;

      if (a_len == 0 || b_len == 0) break;

      // printf("a_seed_begin %d, b_seed_begin %d\n", a_seed_begin, b_seed_begin);

      const int klen = restrictedSize + 2 + 2;
      sType* k1 = &K1[workerId * (klen)];
      sType* k2 = &K2[workerId * (klen)];

      int* _k1 = &k1[2];
      int* _k2 = &k2[2];

      sType partscore = 0;
      // auto tt = TileTelemetry();
      memset(k1, 0, (klen) * sizeof(sType));
      memset(k2, 0, (klen) * sizeof(sType));
      if (a_seed_begin == -1) {
       partscore = 0; 
      } else {
        if (isLeft) {
          partscore = ipumacore::xdrop::xdrop_smart_restricted_extend_left<X, GAP_PENALTY, poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>>&, sType>(
              a, a_len, a_seed_begin,
              b, b_len, b_seed_begin,
              seedLength,
              simMatrix, _k1, _k2, klen);
        } else {
          // Right hand side
          partscore = ipumacore::xdrop::xdrop_smart_restricted_extend_right<X, GAP_PENALTY, poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>>&, sType>(
              a, a_len, a_seed_begin,
              b, b_len, b_seed_begin,
              seedLength,
              simMatrix, _k1, _k2, klen);
        }
      }
      score[n] += partscore;
      // if (isLeft) {
      //   ARange[n] = tt.cycle_counter_diff();
      // } else {
      //   BRange[n] = tt.cycle_counter_diff();
      // }

      // This operation can cause raise conditions, but that will not produce
      // incurrect results. However, double work can be done, which we accept.
      myN = *globalN;
      (*globalN) += 1;
    }
    // printf("AAA: %d\n", (unsigned) tt.cycle_counter_diff());

    // if (workerId == 0)
    //   for (int i = 0; i < 5; i++) {
    //     printf("Score[%d] = %d\n", i, score[i]);
    //   }
    return true;
  }
};

template class SeedExtendRestrictedXDrop<10>;
template class SeedExtendRestrictedXDrop<15>;
template class SeedExtendRestrictedXDrop<20>;
template class SeedExtendRestrictedXDrop<50>;