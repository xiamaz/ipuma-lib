#include <poplar/Vertex.hpp>
#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/AvailableVTypes.h>
#include "poplar/TileConstants.hpp"
#include <print.h>
#include <type_traits>

int max(int a, int b) {
    return a > b ? a : b;
}

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

class XDrop : public poplar::Vertex {
private:
   poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> K1;
   poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> K2;
   poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> K3;
public:
    // Fields
    poplar::Vector<poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>>> simMatrix;
    poplar::Input<size_t> maxNPerTile;
    poplar::Input<int> gapInit;
    poplar::Input<int> gapExt;
    poplar::Input<int> bufSize;
    poplar::Input<int> maxAB;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Seqs;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Meta;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> score;

    bool compute() {

        int gI = *gapInit;
        int gE = *gapExt;

        uint8_t* cSeqs = (uint8_t*) &(Seqs[0]);
        
        for (int n = 0; n < maxNPerTile; ++n) {
            int lastNoGap, prevNoGap;
            int s = 0;
            uint16_t Astart = 0;
            uint16_t Bstart = 0;
            uint16_t Aend = 0;
            uint16_t Bend = 0;

            auto a_len = Meta[4 * n];
            int j_offset = Meta[4 * n + 1];
            auto b_len = Meta[4 * n + 2];
            int i_offset = Meta[4 * n + 3];

            uint8_t* a = cSeqs + j_offset;
            uint8_t* b = cSeqs + i_offset;

            if (a_len == 0 || b_len == 0) break;

            // Algo begin
            int T_prime = 0, T = 0, L = 0, U = 0;

            int M = a_len;
            int N = b_len;

            memset(&(K1[0]), 0, (M + 2) * sizeof(int));
            memset(&(K2[0]), 0, (M + 2) * sizeof(int));
            memset(&(K3[0]), 0, (M + 2) * sizeof(int));

            int* k1 = &K1[1];
            int* k2 = &K2[1];
            int* k3 = &K3[1];

            auto cell_update = [&](int i, int j, int* k1, int* k2, int* k3, int z) {
               auto [index, score] = maxtuple({k2[z] - gE,
                                               k2[z - 1] - gE,
                                               k1[z - 1] + simMatrix[a[i]][b[j]]});
               if (score < T - X) {
                 score = neginf;
               }
               k3[z] = score;
            //    H(i+1, j+1) = score; // DEBUG
               return score;
             };

             auto rotate = [&]() {
               int* t;  // 1->3, 2->1, 3->2 == ((buf + len) % (3*buf))
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
                // C(i + 1, j + 1) = 1; // DEBUG
                T_prime = max(T_prime, score);
              }

              int minL = 9999;
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
            score[n] = T;
        }
        return true;
    }
};