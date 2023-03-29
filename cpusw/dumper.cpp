#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "ipuma.h"
#include "cmd_arguments.hpp"

using json = nlohmann::json;

struct DumperConfig {
	ipu::SWConfig swconfig;
	ipu::LoaderConfig loaderConfig;
	std::string outputSequences;
	std::string outputCmps;
	std::string outputLogan;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DumperConfig, loaderConfig, outputSequences, outputCmps, outputLogan)

void dumpLogan(const ipu::RawSequences& seqs, const ipu::MultiComparisons& cmps, const DumperConfig& config) {
	std::ofstream file(config.outputLogan);
	PLOGI << "Dumping sequences in LOGAN format to " << config.outputLogan;
	for (const auto& mc : cmps) {
		for (const auto& c : mc.comparisons) {
			for (const auto& seed : c.seeds) {
				if (seed.seedAStartPos >= 0 && seed.seedBStartPos >= 0) {
					// std::stringstream ss;
					file << seqs[c.indexA] << "\t" << seed.seedAStartPos << "\t" << seqs[c.indexB] << "\t" << seed.seedBStartPos << "\tn\n";
				}
			}
		}
	}
}

void dumpComparisons(const ipu::RawSequences& seqs, const ipu::MultiComparisons& cmps, const DumperConfig& config) {
	PLOGI << "Dumping sequences to " << config.outputSequences;
	std::ofstream file(config.outputSequences);
	file << json(seqs);

	PLOGI << "Dumping comparisons to " << config.outputCmps;
	std::ofstream cfile(config.outputCmps);
	cfile << json(cmps);
}

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

	cxxopts::Options options("cpusw", "CPU Xdrop implementation");

	options.add_options()
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;

	json configJson = DumperConfig();
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

	DumperConfig config = configJson.get<DumperConfig>();

	PLOGI << "DUMPERCONFIG" << json{config}.dump();
	auto seqdb = config.loaderConfig.getMultiSequences(config.swconfig);
	auto [seqs, mcmps] = seqdb.get();
	PLOGI << "DSSTATS" << ipu::getDatasetStats(seqs, mcmps, config.swconfig.seedLength).dump();

	if (config.outputLogan.size() > 0) {
		dumpLogan(seqs, mcmps, config);
	} else if (config.outputCmps.size() > 0 && config.outputSequences.size() > 0) {
		dumpComparisons(seqs, mcmps, config);
	}

	return 0;

}