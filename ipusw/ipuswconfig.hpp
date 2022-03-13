#ifndef IPU_SW_CONFIG_HPP
#define IPU_SW_CONFIG_HPP

#include "swatlib/swatlib.h"
#include "driver.hpp"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct IpuSwConfig {
	ipu::SWConfig swconfig;
	ipu::IPUAlgoConfig ipuconfig;
	int numDevices = 1; // number of IPU devices to be used
	int numThreads = 1; // number of threads for submitting work to the IPU devices
	bool duplicateDatasets = false;

	IpuSwConfig() : swconfig(SW_CONFIGURATION), ipuconfig(ALGO_CONFIGURATION) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu) : swconfig(sw), ipuconfig(ipu) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu, int numDevices, int numThreads) : swconfig(sw), ipuconfig(ipu), numDevices(numDevices), numThreads(numThreads) {}
};

void to_json(json& j, const IpuSwConfig& c) {
	j = json{
		{"sw", c.swconfig},
		{"ipu", c.ipuconfig},
		{"numDevices", c.numDevices},
		{"numThreads", c.numThreads},
		{"duplicateDatasets", c.duplicateDatasets}
	};
}

void from_json(const json& j, IpuSwConfig& c) {
	j.at("sw").get_to(c.swconfig);
	j.at("ipu").get_to(c.ipuconfig);
	j.at("numDevices").get_to(c.numDevices);
	j.at("numThreads").get_to(c.numThreads);
	j.at("duplicateDatasets").get_to(c.duplicateDatasets);
}
#endif