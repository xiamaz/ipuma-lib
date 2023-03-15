#pragma once
#include "types.h"
#include "ipu_config.h"
#include "sequence_database.h"

namespace ipu {

struct GeneratorConfig {
	int sequenceLength = 0;
	int sequenceCount = 0;
	float similarity = 1.0;
	int seed = 42;
};

class SequenceGenerator : public SequenceDatabase {
	GeneratorConfig config;

	std::vector<std::string> sequences;
	ipu::RawSequences seqs;
	ipu::Comparisons cmps;
public:
	SequenceGenerator(GeneratorConfig, SWConfig);
	virtual std::tuple<RawSequences, ipu::Comparisons> get() override;
};


void to_json(json& j, const GeneratorConfig& c);
void from_json(const json& j, GeneratorConfig& c);
}