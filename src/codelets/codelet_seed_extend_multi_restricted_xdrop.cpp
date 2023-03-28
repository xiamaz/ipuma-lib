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

// const int neginf = -9999;

// typedef int sType;
const int GAP_PENALTY = 1;
typedef float sType;

// typedef int sType;
/*
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
    int64_t now = (((int64_t)now_u) << 32) | now_l;
    int64_t old = (((int64_t)cycle_counter_u) << 32) | cycle_counter_l;
    return (unsigned)now - old;
  }

  void reset() {
    cycle_counter_l = __builtin_ipu_get_scount_l();
    cycle_counter_u = __builtin_ipu_get_scount_u();
  }
};
*/

constexpr unsigned numberOfBits(unsigned x) {
  return x < 2 ? x : 1 + numberOfBits(x >> 1);
}

void inline setZeroPart(sType* ptr, int N) {
  sType negXinf = -9999;
  const float2 zzero = {negXinf, negXinf}; 
  const rptsize_t loopCount = N / 2;
  float2 * o = reinterpret_cast<float2 *>(ptr);
  for (unsigned i = 0; i < loopCount; i++) {
    ipu::store_postinc(&o, zzero, 1);
  }
  ptr[N-1] = negXinf;
}

template <int X>
class SeedExtendRestrictedXDrop : public poplar::MultiVertex {
 private:
  poplar::Output<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR, 8>> K1;
  poplar::Output<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR, 8>> K2;

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
    volatile int* globalWB = &currentN;
    *globalN = 6;

    const ipu::XDropMeta* metas = reinterpret_cast<const ipu::XDropMeta*>(&Meta[0]);
    if (workerId == 0) {
      memset((char*)&score[0], 0, maxNPerTile * sizeof(score[0]));
      // memset((char*)&ARange[0], 0, maxNPerTile * sizeof(ARange[0]));
    }
    for (; myN < maxNPerTile * 2;) {
      // const bool isLeft = myN % 2 == 0;
      const int n = myN / 2;

      const ipu::XDropMeta meta = metas[n];
      const auto a_len = meta.sizeA;
      const auto j_offset = meta.offsetA;
      const auto b_len = meta.sizeB;
      const auto i_offset = meta.offsetB;

      const int32_t a_seed_begin = meta.seedAStartPos;
      const int32_t b_seed_begin = meta.seedBStartPos;

      const uint8_t* a = cSeqs + j_offset;
      const uint8_t* b = cSeqs + i_offset;

      if (a_len == 0 || b_len == 0) break;

      // printf("a_seed_begin %d, b_seed_begin %d\n", a_seed_begin, b_seed_begin);
      const int LR_offsert = 20;

      const int klen = restrictedSize + 2 * LR_offsert;
      sType* k1 = &K1[workerId * (klen)];
      sType* k2 = &K2[workerId * (klen)];

      sType* _k1 = &k1[LR_offsert];
      sType* _k2 = &k2[LR_offsert];

      // comps[myN] = 1;

      // auto tt = TileTelemetry();
      if (a_seed_begin == -1) {
        // if (myN % 2 == 0) {
        //   ARange[n] += tt.cycle_counter_diff();
        // } else {
        //   BRange[n] += tt.cycle_counter_diff();
        // }
        for (int i = 0; i < workerId * 5; i++) {
          asm volatile("nop" ::: "memory");
        }
      } else {
        if (myN % 2 == 0) {
          // memset(k1, 0, (klen) * sizeof(sType));
          // memset(k2, 0, (klen) * sizeof(sType));
          setZeroPart(k1, klen);
          // for (int i = 0; i < klen; i++) {
          //   if (k1[i] >= 0) {
          //     printf("FAILURE %d\n", i);
          //   }
          // }
          setZeroPart(k2, klen);
          // for (int i = 0; i < klen; i++) {
          //   if (k2[i] >= 0) {
          //     printf("FAILURE %d\n", i);
          //   }
          // }
          sType lpartscore = ipumacore::xdrop::xdrop_smart_restricted_extend_left<X, GAP_PENALTY, poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>>&, sType>(
              a, a_len, a_seed_begin,
              b, b_len, b_seed_begin,
              seedLength,
              simMatrix, _k1, _k2, restrictedSize);
          ARange[n] = lpartscore;

        } else {
          // memset(k1, 0, (klen) * sizeof(sType));
          // memset(k2, 0, (klen) * sizeof(sType));
          setZeroPart(k1, klen);
          // for (int i = 0; i < klen; i++) {
          //   if (k1[i] >= 0) {
          //     printf("FAILURE %d\n", i);
          //   }
          // }
          setZeroPart(k2, klen);
          // for (int i = 0; i < klen; i++) {
          //   if (k2[i] >= 0) {
          //     printf("FAILURE %d\n", i);
          //   }
          // }
          sType rpartscore = ipumacore::xdrop::xdrop_smart_restricted_extend_right<X, GAP_PENALTY, poplar::Vector<poplar::Input<poplar::Vector<sType, poplar::VectorLayout::ONE_PTR>>>&, sType>(
                                 a, a_len, a_seed_begin,
                                 b, b_len, b_seed_begin,
                                 seedLength,
                                 simMatrix, _k1, _k2, restrictedSize);
          BRange[n] = rpartscore;
        }

        // Right hand side

        // score[n] = rpartscore + lpartscore;
        // printf("XXXXXXXXXX: %d\n", partscore);
        // partscore = lpartscore + rpartscore + 17;
        // printf("[%d] %d => (%d|%d) :: %d => (%d|%d) :: (%d|%d) = %d\n",n, a_len, a_seed_begin, a_len - a_seed_begin -17, b_len, b_seed_begin, b_len - b_seed_begin - 17 ,lpartscore, rpartscore,lpartscore + rpartscore);
        // printf("%d\n", partscore);
      }
      // printf("RAN\n")
      // score[n] = partscore;
      // if (isLeft) {
      //   ARange[n/2] += tt.cycle_counter_diff();
      // } else {
      //   BRange[n/2] += tt.cycle_counter_diff();
      // }

      // This operation can cause raise conditions, but that will not produce
      // incurrect results. However, double work can be done, which we accept.
      {
        myN = *globalN;
        (*globalWB) = (myN + 1);
      }
    }
    // ARange[workerId] = tt.cycle_counter_diff();
    // printf("AAA: %d\n", (unsigned) tt.cycle_counter_diff());
    // printf("AAA: %d\n", (unsigned) count);

    // printf("AAA: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", comps[0], comps[1], comps[2], comps[3], comps[4], comps[5], comps[6], comps[7], comps[8], comps[9], comps[10], comps[11], comps[12], comps[13], comps[14], comps[15], comps[16], comps[17], comps[18], comps[19]);

    // if (workerId == 0)
    //   for (int i = 0; i < 5; i++) {
    //     printf("Score[%d] = %d\n", i, score[i]);
    //   }
    return true;
  }
};

template class SeedExtendRestrictedXDrop<5>;
template class SeedExtendRestrictedXDrop<10>;
template class SeedExtendRestrictedXDrop<15>;
template class SeedExtendRestrictedXDrop<20>;
template class SeedExtendRestrictedXDrop<50>;
template class SeedExtendRestrictedXDrop<100>;
template class SeedExtendRestrictedXDrop<500>;
// template class SeedExtendRestrictedXDrop<20>;
// template class SeedExtendRestrictedXDrop<25>;
