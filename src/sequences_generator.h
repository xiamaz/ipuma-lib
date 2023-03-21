#pragma once
#include "types.h"
#include "ipu_config.h"
#include "sequence_database.h"

namespace ipu {

struct GeneratorConfig {
	int generatorCount = 0;
	int generatorSeqLen = 0;
	float generatorSimilarity = 1.0;
	int generatorSeed = 42;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(GeneratorConfig, generatorCount, generatorSeqLen, generatorSimilarity, generatorSeed);

template<typename C>
SequenceDatabase<C> loadSequences(GeneratorConfig& config, SWConfig& swconfig);

}