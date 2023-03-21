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

size_t getNumberComparisons(const ipu::MultiComparisons& mcmps) {
	size_t totalCmps = 0;
	for (const auto& m : mcmps) {
		totalCmps += m.comparisons.size();
	}
	return totalCmps;
}

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

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

	ipu::SequenceDatabase<ipu::MultiComparison> seqdb = config.loaderconfig.getMultiSequences(config.swconfig);
	auto [seqs, mcmps] = seqdb.get();
	PLOGI << ipu::getDatasetStats(seqs, mcmps).dump();

	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, config.numDevices);

	std::vector<ipu::Batch> batches;
	if (config.decomposeMulticomparisons) {
		PLOGI << "Decompose multicomparisons to individual comparisons";
		auto cmps = ipu::convertToComparisons(mcmps);
		batches = ipu::create_batches(seqs, cmps, config.ipuconfig, config.swconfig);
	} else {
		PLOGI << "Use multicomparisons for creating batches";
		batches = ipu::create_batches(seqs, mcmps, config.ipuconfig, config.swconfig);
	}

	std::vector<int32_t> scores(getNumberComparisons(mcmps));
	swatlib::TickTock time;
	double gcells = 0;
	time.tick();
	std::vector<ipu::batchaffine::Job*> jobs(batches.size());
	for (int i = 0; i < batches.size(); ++i) {
    jobs[i] = driver.async_submit(&batches[i]);
	}
	for (int i = 0; i < jobs.size(); ++i) {
		driver.blocking_join(*jobs[i]);
		PLOGI << "Received batch " << i+1 << " / " << jobs.size();
		gcells += batches[i].cellCount / 1e9;
	}
	time.tock();
	for (int b = 0; b < batches.size(); ++b) {
    auto result = batches[b].get_result();
    // delete jobs[i];
		for (int i = 0; i < batches[b].origin_comparison_index.size(); ++i) {
			auto [orig_i, orig_seed] = ipu::unpackOriginIndex(batches[b].origin_comparison_index[i]);
			if (orig_i >= 0) {
				int lpartScore = result.a_range_result[i];
				int rpartScore = result.b_range_result[i];
				// PLOGF << orig_i << ":" << batch.origin_comparison_index[i] << " " << lpartScore << " " << rpartScore;
				scores[orig_i] = std::max(lpartScore + rpartScore + config.swconfig.seedLength, scores[orig_i]);
			}
		}
	}

	auto totalTimeMs = time.duration();
	auto gcups = gcells / totalTimeMs * 1e3;

	PLOGI << json{
		{"time_ms", totalTimeMs},
		{"gcups", gcups},
		{"gigacells", gcells},
		{"comparisons", mcmps.size()},
	}.dump();

	// PLOGE << json{scores}.dump();

	if (config.output != "") {
		std::ofstream ofile(config.output);
		ofile << json{scores};
	}
  return 0;
}