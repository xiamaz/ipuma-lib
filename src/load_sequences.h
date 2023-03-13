#pragma once

#include<tuple>
#include<istream>
#include <nlohmann/json.hpp>

#include <plog/Log.h>

#include "types.h"

using json = nlohmann::json;

namespace ipu {
json getDatasetStats(const ipu::RawSequences& seqs, const ipu::Comparisons& cmps);

class SequenceData {
	std::vector<std::string> sequences;
	ipu::RawSequences seqs;
	ipu::Comparisons cmps;

public:
	SequenceData(std::string, std::string, std::string, std::string);
	// SequenceData(std::string, std::string, std::string, std::string, std::string, std::string);

	std::tuple<ipu::RawSequences, ipu::Comparisons> get();
};
}