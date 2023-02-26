
#ifndef __DRIVER_H__
#define __DRIVER_H__


#ifndef KLIGN_IPUS_LOCAL
#define KLIGN_IPUS_LOCAL 1
#endif

#include "ipu_batch_affine.h"

ipu::batchaffine::SWAlgorithm* getDriver();
void init_single_ipu(ipu::SWConfig config, ipu::IPUAlgoConfig algoconfig);



#endif // __DRIVER_H__