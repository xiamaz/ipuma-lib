#ifndef IPU_SW_CONFIG_HPP
#define IPU_SW_CONFIG_HPP

#include "swatlib/swatlib.h"
#include "driver.hpp"
#include "sequences_generator.h"
#include "sequences_loader.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct IpuSwConfig {
	ipu::SWConfig swconfig;
	ipu::IPUAlgoConfig ipuconfig;
	ipu::LoaderConfig loaderconfig;
	int numDevices = 1; // number of IPU devices to be used
	int numThreads = 1; // number of threads for submitting work to the IPU devices
	std::string output = "";
	bool decomposeMulticomparisons = false;

	IpuSwConfig() {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu) : swconfig(sw), ipuconfig(ipu) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu, int numDevices, int numThreads) : swconfig(sw), ipuconfig(ipu), numDevices(numDevices), numThreads(numThreads) {}
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(IpuSwConfig, swconfig, ipuconfig, loaderconfig, numDevices, numThreads, output, decomposeMulticomparisons);

#endif