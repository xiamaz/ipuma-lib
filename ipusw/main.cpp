#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "ipuswconfig.hpp"
#include "ipuma.h"
#include "cmd_arguments.hpp"

using json = nlohmann::json;

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::verbose, &consoleAppender);

	cxxopts::Options options("ipusw", "IPU Smith Waterman Binary");

	options.add_options()
		("hSequencePath", "Sequences File (either fa or txt)", cxxopts::value<std::string>())
		("vSequencePath", "Sequences File (either fa or txt)", cxxopts::value<std::string>())
		("hSeedPath", "Seed Position File (txt)", cxxopts::value<std::string>())
		("vSeedPath", "Seed Position File (txt)", cxxopts::value<std::string>())
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;

	json configJson = IpuSwConfig();
	addArguments(configJson, options, "");

	options.positional_help("[hSequences] [vSequences] [hSeed] [vSeed]");
	options.parse_positional({"hSequencePath", "vSequencePath", "hSeedPath", "vSeedPath"});

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

	IpuSwConfig config = configJson.get<IpuSwConfig>();

	std::string hPath = result["hSequencePath"].as<std::string>();
	std::string vPath = result["vSequencePath"].as<std::string>();
	std::string hSeedPath = result["hSeedPath"].as<std::string>();
	std::string vSeedPath = result["vSeedPath"].as<std::string>();

	PLOGI << "IPUSWCONFIG" << json{config}.dump();

	// run_comparison(config, refPath, queryPath);
	auto seqdb = ipu::SequenceData(hPath, vPath, hSeedPath, vSeedPath);
	auto [seqs, cmps] = seqdb.get();
	PLOGI << ipu::getDatasetStats(seqs, cmps).dump();

	ipu::MultiComparisons mcmps;
	for (const auto& cmp : cmps) {
		mcmps.push_back({{cmp}, 17});
	}

	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, config.numDevices);

	auto batches = ipu::create_batches(seqs, mcmps, config.ipuconfig, config.swconfig);


  std::vector<ipu::BlockAlignmentResults> results;
  for (auto& batch : batches) {
    ipu::batchaffine::Job* j = driver.async_submit(&batch);
    assert(batch.cellCount > 0);
    assert(batch.dataCount > 0);
    driver.blocking_join(*j);
    results.push_back(batch.get_result());
    delete j;
  }
  return 0;
}