#include "swatlib/swatlib.h"
#include "ipuma.h"
#include "driver.hpp"
#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

using json = nlohmann::json;


struct IpuSwConfig {
	ipu::SWConfig swconfig;
	ipu::IPUAlgoConfig ipuconfig;

	IpuSwConfig() : swconfig(SW_CONFIGURATION), ipuconfig(ALGO_CONFIGURATION) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu) : swconfig(sw), ipuconfig(ipu) {}

	IpuSwConfig(json data) {
		swconfig = {
		};
		ipuconfig = {
		};
	}
};

int main(int argc, char** argv) {
	cxxopts::Options options("ipusw", "IPU Smith Waterman Binary");

	options.add_options()
		("reference", "Reference File (either fa or txt)", cxxopts::value<std::string>())
		("query", "Query File (either fa or txt)", cxxopts::value<std::string>())
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;
	options.positional_help("[reference_file] [query_file]");
	options.parse_positional({"reference", "query", ""});

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << "\n";
		exit(0);
	}
}