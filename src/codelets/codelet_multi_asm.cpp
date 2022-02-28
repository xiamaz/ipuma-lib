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
class MultiSWAffineAsm : public poplar::MultiVertex {
private:
   poplar::Output<poplar::Vector<float, poplar::VectorLayout::ONE_PTR, 8>> C;
   poplar::Output<poplar::Vector<float, poplar::VectorLayout::ONE_PTR, 8>> bG;
public:
    // Fields
    poplar::Vector<poplar::Input<poplar::Vector<float, poplar::VectorLayout::ONE_PTR>>, poplar::VectorLayout::ONE_PTR> simMatrix;
    poplar::Input<int> maxNPerTile;
    poplar::Input<float> gapInit;
    poplar::Input<float> gapExt;
    poplar::Input<int> bufSize;
    poplar::Input<int> maxAB;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> Seqs;
    poplar::Input<poplar::Vector<int, poplar::VectorLayout::ONE_PTR, 8>> Meta;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> score;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> ARange;
    poplar::Output<poplar::Vector<int, poplar::VectorLayout::ONE_PTR>> BRange;
    poplar::Input<bool> forwardOnly;

                IS_EXTERNAL_CODELET(1);

    bool compute(unsigned workerId) {
        return true;
    }
};