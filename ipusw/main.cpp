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
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;

	json configJson = IpuSwConfig();
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

	IpuSwConfig config = configJson.get<IpuSwConfig>();

	PLOGI << "IPUSWCONFIG" << json{config}.dump();

	auto seqdb = config.getSequences();
	auto [seqs, cmps] = seqdb->get();
	PLOGI << ipu::getDatasetStats(seqs, cmps).dump();
	if (cmps.size() == 0) {
		PLOGF << "No comparisons defined";
		exit(1);
	}

	ipu::MultiComparisons mcmps;
	for (const auto& cmp : cmps) {
		mcmps.push_back({{cmp}, config.swconfig.seedLength});
	}

	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, config.numDevices);

	auto batches = ipu::create_batches(seqs, cmps, config.ipuconfig, config.swconfig);


	std::vector<int32_t> scores(cmps.size());
  for (auto& batch : batches) {
    ipu::batchaffine::Job* j = driver.async_submit(&batch);
    assert(batch.cellCount > 0);
    assert(batch.dataCount > 0);
    driver.blocking_join(*j);
    auto result = batch.get_result();
    delete j;

		PLOGE << json{result.a_range_result}.dump();

		// results parsing
		for (int i = 0; i < batch.origin_comparison_index.size(); ++i) {
			auto [orig_i, orig_seed] = ipu::unpackOriginIndex(batch.origin_comparison_index[i]);
			if (orig_i >= 0) {
				int lpartScore = result.a_range_result[i];
				int rpartScore = result.a_range_result[i];
				PLOGF << orig_i << ":" << batch.origin_comparison_index[i] << " " << lpartScore << " " << rpartScore;
				scores[orig_i] = std::max(lpartScore + rpartScore, scores[orig_i]);
			}
		}
  }

	PLOGE << json{scores}.dump();

	if (config.output != "") {
		std::ofstream ofile(config.output);
		ofile << json{scores};
	}
  return 0;
}