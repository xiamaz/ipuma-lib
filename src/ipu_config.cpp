#include "ipu_config.h"

namespace ipu {
int IPUAlgoConfig::getBufsize32b() const { return std::ceil(static_cast<double>(bufsize) / 4.0); }

int IPUAlgoConfig::getTotalNumberOfComparisons() const { return tilesUsed * maxBatches; }

int IPUAlgoConfig::getMetaBufferSize32b() const { return getTotalNumberOfComparisons() * 4; };

int IPUAlgoConfig::getTotalBufsize32b() const { return tilesUsed * getBufsize32b(); }

int IPUAlgoConfig::getInputBufferSize32b() const { return getTotalBufsize32b() + getMetaBufferSize32b(); }
}