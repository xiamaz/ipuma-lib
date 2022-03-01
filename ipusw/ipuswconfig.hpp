#ifndef IPU_SW_CONFIG_HPP
#define IPU_SW_CONFIG_HPP

#include "swatlib/swatlib.h"
#include "driver.hpp"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct IpuSwConfig {
	ipu::SWConfig swconfig;
	ipu::IPUAlgoConfig ipuconfig;
	int numDevices = 1;

	IpuSwConfig() : swconfig(SW_CONFIGURATION), ipuconfig(ALGO_CONFIGURATION) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu) : swconfig(sw), ipuconfig(ipu) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu, int numDevices) : swconfig(sw), ipuconfig(ipu), numDevices(numDevices) {}
};

void to_json(json& j, const IpuSwConfig& c) {
	j = json{
		{"sw", c.swconfig},
		{"ipu", c.ipuconfig},
		{"numDevices", c.numDevices}
	};
}

void from_json(const json& j, IpuSwConfig& c) {
	j.at("sw").get_to(c.swconfig);
	j.at("ipu").get_to(c.ipuconfig);
	j.at("numDevices").get_to(c.numDevices);
}
#endif