#ifndef IPU_SW_CONFIG_HPP
#define IPU_SW_CONFIG_HPP

#include "swatlib/swatlib.h"
#include "driver.hpp"
#include "sequences_generator.h"
#include "load_sequences.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct IpuSwConfig {
	ipu::SWConfig swconfig;
	ipu::IPUAlgoConfig ipuconfig;
	ipu::GeneratorConfig generatorconfig;
	ipu::SequenceConfig sequenceconfig;
	bool generateSequences = false;
	int numDevices = 1; // number of IPU devices to be used
	int numThreads = 1; // number of threads for submitting work to the IPU devices
	std::string output = "";
	int repeat = 1; // number of repeated runs

	IpuSwConfig() {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu) : swconfig(sw), ipuconfig(ipu) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu, int numDevices, int numThreads) : swconfig(sw), ipuconfig(ipu), numDevices(numDevices), numThreads(numThreads) {}

	std::unique_ptr<ipu::SequenceDatabase> getSequences() {
		if (generateSequences) {
			return std::unique_ptr<ipu::SequenceDatabase>(new ipu::SequenceGenerator(generatorconfig, swconfig));
		} else {
			return std::unique_ptr<ipu::SequenceDatabase>(new ipu::SequenceData(sequenceconfig));
		}
	}
};

void to_json(json& j, const IpuSwConfig& c) {
	j = json{
		{"sw", c.swconfig},
		{"ipu", c.ipuconfig},
		{"generator", c.generatorconfig},
		{"sequence", c.sequenceconfig},
		{"generateSequences", c.generateSequences},
		{"numDevices", c.numDevices},
		{"numThreads", c.numThreads},
		{"output", c.output},
		{"repeat", c.repeat}
	};
}

void from_json(const json& j, IpuSwConfig& c) {
	j.at("sw").get_to(c.swconfig);
	j.at("ipu").get_to(c.ipuconfig);
	j.at("numDevices").get_to(c.numDevices);
	j.at("generator").get_to(c.generatorconfig);
	j.at("sequence").get_to(c.sequenceconfig);
	j.at("generateSequences").get_to(c.generateSequences);
	j.at("numThreads").get_to(c.numThreads);
	j.at("output").get_to(c.output);
	j.at("repeat").get_to(c.repeat);
}
#endif