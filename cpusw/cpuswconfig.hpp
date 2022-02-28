#ifndef CPU_SW_CONFIG_HPP
#define CPU_SW_CONFIG_HPP

#include <nlohmann/json.hpp>
#include "swatlib/swatlib.h"
#include "ipu_config.h"

using json = nlohmann::json;

namespace cpu {
	static const std::vector<std::string> algoNames = {"stripedsw", "seqan3"};
	enum class Algo {stripedsw, seqan3};
	void to_json(json& j, const Algo& a) {
		j = algoNames.at(static_cast<int>(a));
	}
	void from_json(const json& j, Algo& a) {
		const std::string& s = j;
		for (int i = 0; i < algoNames.size(); ++i) {
			if (algoNames[i] == s) { a = static_cast<Algo>(i); return; }
		}
		throw std::runtime_error("Invalid cpu algo type string: " + s);
	}

	struct AlgoConfig {
		int threads = 1;
		Algo algo = Algo::stripedsw;
	};

	void to_json(json& j, const AlgoConfig& c) {
		j = {{"threads", c.threads}, {"algo", c.algo}};
	}

	void from_json(const json& j, AlgoConfig& c) {
		j.at("threads").get_to(c.threads);
		j.at("algo").get_to(c.algo);
	}
}

struct CpuSwConfig {
	ipu::SWConfig swconfig;
	cpu::AlgoConfig algoconfig;
};

void to_json(json& j, const CpuSwConfig& c) {
	j = json{
		{"sw", c.swconfig},
		{"cpu", c.algoconfig}
	};
}

void from_json(const json& j, CpuSwConfig& c) {
	j.at("sw").get_to(c.swconfig);
	j.at("cpu").get_to(c.algoconfig);
}

#endif