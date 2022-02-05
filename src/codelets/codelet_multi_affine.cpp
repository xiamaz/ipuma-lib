#include <poplar/Vertex.hpp>
#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/AvailableVTypes.h>
#include "poplar/TileConstants.hpp"
#include <print.h>
#include <type_traits>

inline int max(int a, int b) {
    return a > b ? a : b;
}

/**
 * Single SW operations for cell i, j
 * based on Fig. 1, Wozniak et al 1997
 */
class MultiSWAffine : public poplar::MultiVertex {
private:
    // poplar::Vector<poplar::Vector<int, poplar::VectorLayout::ONE_PTR, 8>, poplar::VectorLayout::ONE_PTR> C;
    // poplar::Vector<poplar::Vector<int, poplar::VectorLayout::ONE_PTR, 8>, poplar::VectorLayout::ONE_PTR> bG;
    poplar::Vector<int, poplar::VectorLayout::ONE_PTR, 8> C;
    poplar::Vector<int, poplar::VectorLayout::ONE_PTR, 8> bG;
public:
    // Fields
    poplar::Vector<poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>>, poplar::VectorLayout::ONE_PTR> simMatrix;
    poplar::Input<int> maxNPerTile;
    poplar::Input<int> gapInit;
    poplar::Input<int> gapExt;
    poplar::Input<int> bufSize;
    poplar::Input<int> maxAB;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> A;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> B;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Alen;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Blen;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> score;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> mismatches;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> ARange;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> BRange;

    bool compute(unsigned workerId) {
        int gI = *gapInit;
        int gE = *gapExt;

        auto* wC = &(C[0]) + workerId * maxAB;
        auto* wbG = &(bG[0]) + workerId * maxAB;

        uint8_t* cA = (uint8_t*) &(A[0]);
        uint8_t* cB = (uint8_t*) &(B[0]);

        for (int n = workerId; n < maxNPerTile; n += MultiVertex::numWorkers()) {
            int lastNoGap, prevNoGap;
            int s = 0;
            uint16_t Astart = 0;
            uint16_t Bstart = 0;
            uint16_t Aend = 0;
            uint16_t Bend = 0;

            auto a_len = Alen[2*n];
            auto b_len = Blen[2*n];
            auto j_offset = Alen[2*n+1];
            auto i_offset = Blen[2*n+1];
            if (a_len == 0 || b_len == 0) break;

            memset(&(wC[0]), 0, maxAB * sizeof(int));
            memset(&(wbG[0]), 0, maxAB * sizeof(int));

            // forward pass
            for (int i = 0; i < b_len; ++i) {
                int aGap;
                lastNoGap = prevNoGap = 0;
                aGap = gapInit;
                for (unsigned j = 0; j < a_len; ++j) {
                    aGap = max(lastNoGap + gI + gE, aGap + gE);
                    wbG[j] = max(wC[j] + gI + gE, wbG[j] + gE);

                    lastNoGap = max(prevNoGap + simMatrix[cA[j_offset + j]][cB[i_offset + i]], aGap);
                    lastNoGap = max(lastNoGap, wbG[j]);
                    lastNoGap = max(lastNoGap, 0);
                    prevNoGap = wC[j];
                    wC[j] = lastNoGap;
                    if (lastNoGap > s) {
                        Aend = j;
                        Bend = i;
                        s = lastNoGap;
                    }
                }
            }
            score[n] = s;

            s = 0;

            memset(&(wC[0]), 0, maxAB * sizeof(int));
            memset(&(wbG[0]), 0, maxAB * sizeof(int));

            // reverse pass
            for (int i = Bend; i >= 0; --i) {
                int aGap;
                lastNoGap = prevNoGap = 0;
                aGap = gapInit;
                for (int j = Aend; j >= 0; --j) {
                    aGap = max(lastNoGap + gI + gE, aGap + gE);
                    wbG[j] = max(wC[j] + gI + gE, wbG[j] + gE);
                    lastNoGap = max(prevNoGap + simMatrix[cA[j_offset + j]][cB[i_offset + i]], aGap);
                    lastNoGap = max(lastNoGap, wbG[j]);
                    lastNoGap = max(lastNoGap, 0);
                    prevNoGap = wC[j];
                    wC[j] = lastNoGap;
                    if (lastNoGap > s) {
                        Astart = j;
                        Bstart = i;
                        s = lastNoGap;
                    }
                }
            }

            uint16_t range[2];
            range[0] = Astart;
            range[1] = Aend;
            ARange[n] = *reinterpret_cast<uint32_t*>(range);
            range[0] = Bstart;
            range[1] = Bend;
            BRange[n] = *reinterpret_cast<uint32_t*>(range);
        }
        return true;
    }
};