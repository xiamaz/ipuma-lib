
#ifndef __DRIVER_H__
#define __DRIVER_H__

#ifndef KLIGN_IPU_MAXAB_SIZE
#define KLIGN_IPU_MAXAB_SIZE 200
#endif

#ifndef KLIGN_IPU_TILES
#define KLIGN_IPU_TILES 8832
#endif

#ifndef KLIGN_IPUS_LOCAL
#define KLIGN_IPUS_LOCAL 1
#endif

#ifndef KLIGN_IPU_MAX_BATCHES
#define KLIGN_IPU_MAX_BATCHES 300
#endif

#ifndef KLIGN_IPU_BUFSIZE
#define KLIGN_IPU_BUFSIZE 30000
#endif

#include "ipu_batch_affine.h"

static const ipu::SWConfig SW_CONFIGURATION = {
              .gapInit = -(ALN_GAP_OPENING_COST - ALN_GAP_EXTENDING_COST),
              .gapExtend = -ALN_GAP_EXTENDING_COST,
              .matchValue = ALN_MATCH_SCORE,
              .mismatchValue = -ALN_MISMATCH_COST,
              .ambiguityValue = -ALN_AMBIGUITY_COST,
              .similarity = swatlib::Similarity::nucleicAcid,
              .datatype = swatlib::DataType::nucleicAcid,
};

static const ipu::batchaffine::IPUAlgoConfig ALGO_CONFIGURATION = {
            KLIGN_IPU_TILES,
            KLIGN_IPU_MAXAB_SIZE,
            KLIGN_IPU_MAX_BATCHES,
            KLIGN_IPU_BUFSIZE,
            ipu::batchaffine::VertexType::cpp,
            ipu::batchaffine::partition::Algorithm::fillFirst
};

ipu::batchaffine::SWAlgorithm* getDriver();
void init_single_ipu(ipu::SWConfig config, ipu::batchaffine::IPUAlgoConfig algoconfig);



#endif // __DRIVER_H__