#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "ipuswconfig.hpp"
#include "ipuma.h"

using json = nlohmann::json;

// recursively parse arguments
void parseArguments(json& dict, const cxxopts::ParseResult& result) {
	for (const auto& [k, v] : dict.items()) {
		switch(v.type()) {
		case json::value_t::object:
			parseArguments(v, result);
			break;
		case json::value_t::number_integer:
			dict[k] = result[k].as<int>();
			break;
		case json::value_t::number_float:
			dict[k] = result[k].as<double>();
			break;
		case json::value_t::string:
			dict[k] = result[k].as<std::string>();
			break;
		case json::value_t::boolean:
			dict[k] = result[k].as<bool>();
			break;
		default:
			throw std::runtime_error("unsupported type");
			break;
		}
	}
}

void addArguments(const json& dict, cxxopts::Options& options, std::string gname) {
	std::string strval;
	for (const auto& [k, v] : dict.items()) {
		switch (v.type()) {
			case json::value_t::number_integer:
				strval = std::to_string(v.get<int>());
				options.add_option(gname, "", k, "", cxxopts::value<int>()->default_value(strval), strval);
				break;
			case json::value_t::number_float:
				strval = std::to_string(v.get<double>());
				options.add_option(gname, "", k, "", cxxopts::value<double>()->default_value(strval), strval);
				break;
			case json::value_t::string:
				options.add_option(gname, "", k, "", cxxopts::value<std::string>()->default_value(v), v);
				break;
			case json::value_t::boolean:
				strval = std::to_string(v.get<bool>());
				options.add_option(gname, "", k, "", cxxopts::value<bool>()->default_value(strval), strval);
				break;
			case json::value_t::object:
				addArguments(v, options, k);
				break;
			default:
				throw std::runtime_error("unsupported type");
				break;
		}
	}
}

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

	cxxopts::Options options("ipusw", "IPU Smith Waterman Binary");

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

	std::string hPath = result["hSequencePath"].as<std::string>();
	std::string vPath = result["vSequencePath"].as<std::string>();
	std::string hSeedPath = result["hSeedPath"].as<std::string>();
	std::string vSeedPath = result["vSeedPath"].as<std::string>();

	PLOGI << "IPUSWCONFIG" << json{config}.dump();

	// run_comparison(config, refPath, queryPath);
	auto [seqs, cmps] = ipu::prepareComparisons(hPath, vPath, hSeedPath, vSeedPath);
	PLOGI << ipu::getDatasetStats(seqs, cmps).dump();

	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, config.numDevices);

	auto batches = driver.create_batches(seqs, cmps);

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