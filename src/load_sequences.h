#pragma once

#include<tuple>
#include<istream>
#include <nlohmann/json.hpp>

#include <plog/Log.h>

#include "types.h"

using json = nlohmann::json;

namespace ipu {
json getDatasetStats(const ipu::RawSequences& seqs, const ipu::Comparisons& cmps);

std::tuple<ipu::RawSequences, ipu::Comparisons> prepareComparisons(std::string seqH, std::string seqV, std::string seedH, std::string seedV);
}