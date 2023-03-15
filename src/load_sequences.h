#pragma once

#include<tuple>
#include<istream>
#include <nlohmann/json.hpp>

#include <plog/Log.h>

#include "types.h"
#include "sequence_database.h"

using json = nlohmann::json;

namespace ipu {
json getDatasetStats(const ipu::RawSequences& seqs, const ipu::Comparisons& cmps);

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

void from_json(const json& j, SequenceConfig& c);

void to_json(json& j, const SequenceConfig& c);

class SequenceData : public SequenceDatabase {
	std::vector<std::string> sequences;
	ipu::RawSequences seqs;
	ipu::Comparisons cmps;

public:
	SequenceData(const SequenceConfig& c);
	// SequenceData(std::string, std::string, std::string, std::string, std::string, std::string);

	virtual std::tuple<ipu::RawSequences, ipu::Comparisons> get() override;
};
}