#include "driver.hpp"
#include "ipu_base.h"
#include "ipu_batch_affine.h"

static ipu::batchaffine::SWAlgorithm *ipu_driver = nullptr;

ipu::batchaffine::SWAlgorithm* getDriver() {
      assert(ipu_driver != nullptr);
        return ipu_driver;
}

void init_single_ipu(ipu::SWConfig config, ipu::IPUAlgoConfig algoconfig) {
  if (ipu_driver == NULL) {
    ipu_driver =
        new ipu::batchaffine::SWAlgorithm(config, algoconfig, true, 10);
  }
}
