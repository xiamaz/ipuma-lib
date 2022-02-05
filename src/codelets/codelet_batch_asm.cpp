#include <poplar/Vertex.hpp>
#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/AvailableVTypes.h>
#include <print.h>
#include "ExternalCodelet.hpp"

/**
 * Single SW operations for cell i, j
 * based on Fig. 1, Wozniak et al 1997
 */
class SWAffineAsm : public poplar::Vertex {
private:
    poplar::Vector<float, poplar::VectorLayout::ONE_PTR, 8> C;
    poplar::Vector<float, poplar::VectorLayout::ONE_PTR, 8> bG;
public:
    // Fields
    poplar::Vector<poplar::Input<poplar::Vector<float, poplar::VectorLayout::ONE_PTR>>, poplar::VectorLayout::ONE_PTR> simMatrix;
    poplar::Input<int> maxNPerTile;
    poplar::Input<float> gapInit;
    poplar::Input<float> gapExt;
    poplar::Input<int> bufSize;
    poplar::Input<int> maxAB;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> A;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> B;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR, 8>> Alen;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR, 8>> Blen;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> score;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> mismatches;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> ARange;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> BRange;

		IS_EXTERNAL_CODELET(1);

    bool compute() {
        return true;
    }
};