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
#include <omp.h>

#define LAZY_GENERATION_THRES 750

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
	PLOGI << ipu::getDatasetStats(seqs, mcmps, config.swconfig.seedLength).dump();

	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, config.numDevices);

	std::vector<ipu::partition::BatchMapping> mappings;

	// std::vector<ipu::Batch> batches;
	if (config.decomposeMulticomparisons) {
		PLOGI << "Decompose multicomparisons to individual comparisons";
		auto cmps = ipu::convertToComparisons(mcmps);
		mappings = ipu::partition::mapBatches(config.ipuconfig, seqs, cmps);
		// batches = ipu::create_batches(seqs, cmps, config.ipuconfig, config.swconfig);
	} else {
		PLOGI << "Use multicomparisons for creating batches";
		mappings = ipu::partition::mapBatches(config.ipuconfig, seqs, mcmps);
		// batches = ipu::create_batches(seqs, mcmps, config.ipuconfig, config.swconfig);
	}

	std::vector<int32_t> scores(getNumberComparisons(mcmps));
	swatlib::TickTock time;
	double gcells = 0;
	std::vector<ipu::batchaffine::Job*> jobs(mappings.size());
	std::vector<ipu::Batch> batches(mappings.size());
	bool lazyGeneration = false;

	if (mappings.size() > 750) {
		lazyGeneration = true;
		PLOGI << "Generate batches lazy, as batches " << mappings.size() << " > " << LAZY_GENERATION_THRES;
	}

	if (!lazyGeneration) {
		#pragma omp parallel for
		for (int i = 0; i < mappings.size(); ++i) {
			batches[i] = ipu::create_batch(mappings[i], seqs, config.ipuconfig, config.swconfig);
		}
	}

	int progress = 0;

	time.tick();
	#pragma omp parallel for
	for (int i = 0; i < mappings.size(); ++i) {
		ipu::Batch batch;
		if (lazyGeneration) {
			batch = ipu::create_batch(mappings[i], seqs, config.ipuconfig, config.swconfig);
		} else {
			batch = std::move(batches[i]);
		}
    jobs[i] = driver.async_submit(&batch);
		driver.blocking_join(*jobs[i]);

		#pragma omp critical
		{
			progress++;
			gcells += batch.cellCount / 1e9;
			PLOGI << "Received batch " << progress << " / " << jobs.size();
		}

    auto result = batch.get_result();
    delete jobs[i];
		for (int ii = 0; ii < batch.origin_comparison_index.size(); ++ii) {
			auto [orig_i, orig_seed] = ipu::unpackOriginIndex(batch.origin_comparison_index[ii]);
			if (orig_i >= 0) {
				int lpartScore = result.a_range_result[ii];
				int rpartScore = result.b_range_result[ii];
				// PLOGF << orig_i << ":" << batch.origin_comparison_index[i] << " " << lpartScore << " " << rpartScore;
				scores[orig_i] = std::max(lpartScore + rpartScore + config.swconfig.seedLength, scores[orig_i]);
			}
		}
	}
	time.tock();

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