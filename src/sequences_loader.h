#pragma once

#include<tuple>
#include<istream>
#include <nlohmann/json.hpp>

#include <plog/Log.h>

#include "types.h"
#include "sequence_database.h"
#include "sequences_generator.h"

using json = nlohmann::json;

namespace ipu {
template<typename C>
json getDatasetStats(const ipu::RawSequences& seqs, const std::vector<C>& cmps, int seedLen);

struct JsonSequenceConfig {
	std::string comparisons;
	std::string sequences;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(JsonSequenceConfig, comparisons, sequences);

struct SequenceConfig {
	std::string seqsH;
	std::string seqsV;
	std::string seedsH1;
	std::string seedsV1;
	std::string seedsH2;
	std::string seedsV2;
	std::string complexity1;
	std::string complexity2;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SequenceConfig, seqsH, seqsV, seedsH1, seedsV1, seedsH2, seedsV2, complexity1, complexity2);

template<typename C>
class SequenceData : public SequenceDatabase<C> {
	std::vector<std::string> sequences;
public:
	SequenceData(const SequenceConfig& c);
	// SequenceData(std::string, std::string, std::string, std::string, std::string, std::string);

	SequenceData(const JsonSequenceConfig& c);
};

template<typename C>
SequenceDatabase<C> loadSequences(SequenceConfig& config, SWConfig& swconfig);

template<typename C>
SequenceDatabase<C> loadSequences(JsonSequenceConfig& config);

struct LoaderConfig {
	bool useMulticomparison = false;
	SequenceConfig seqConfig;
	JsonSequenceConfig jsConfig;
	GeneratorConfig genConfig;

	ipu::SequenceDatabase<Comparison> getSequences(SWConfig swconfig);
	ipu::SequenceDatabase<MultiComparison> getMultiSequences(SWConfig swconfig);
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(LoaderConfig, seqConfig, jsConfig, genConfig, useMulticomparison);

}