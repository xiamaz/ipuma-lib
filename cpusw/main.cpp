#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "cpuswconfig.hpp"
#include "run_comparison.hpp"

#include "ipuma.h"
#include "cmd_arguments.hpp"
#include "alignment_seqan.hpp"
#include "alignment_genometools.hpp"

using json = nlohmann::json;

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

	cxxopts::Options options("cpusw", "CPU Xdrop implementation");

	options.add_options()
		("hSequencePath", "Sequences File (either fa or txt)", cxxopts::value<std::string>())
		("vSequencePath", "Sequences File (either fa or txt)", cxxopts::value<std::string>())
		("hSeedPath", "Seed Position File (txt)", cxxopts::value<std::string>())
		("vSeedPath", "Seed Position File (txt)", cxxopts::value<std::string>())
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;

	options.positional_help("[hSequences] [vSequences] [hSeed] [vSeed]");
	options.parse_positional({"hSequencePath", "vSequencePath", "hSeedPath", "vSeedPath"});

	json configJson = CpuSwConfig();
	addArguments(configJson, options, "");

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << "\n";
		exit(0);
	}

	if (result.count("config")) {
		std::string configPath = result["config"].as<std::string>();
		std::ifstream cf(configPath);
		json cj;
		cf >> cj;
		configJson = cj;
	}

	parseArguments(configJson, result);

	CpuSwConfig config = configJson.get<CpuSwConfig>();

	std::string hPath = result["hSequencePath"].as<std::string>();
	std::string vPath = result["vSequencePath"].as<std::string>();
	std::string hSeedPath = result["hSeedPath"].as<std::string>();
	std::string vSeedPath = result["vSeedPath"].as<std::string>();

	PLOGI << "CPUSWCONFIG" << json{config}.dump();

	test();
	return 0;

	auto [seqs, cmps] = ipu::prepareComparisons(hPath, vPath, hSeedPath, vSeedPath);
	PLOGI << ipu::getDatasetStats(seqs, cmps).dump();

	runAlignment(seqs, cmps);
}