#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "cpuswconfig.hpp"
#include "run_comparison.hpp"

using json = nlohmann::json;

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

	cxxopts::Options options("ipusw", "IPU Smith Waterman Binary");

	options.add_options()
		("reference", "Reference File (either fa or txt)", cxxopts::value<std::string>())
		("query", "Query File (either fa or txt)", cxxopts::value<std::string>())
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;

	options.positional_help("[reference_file] [query_file]");
	options.parse_positional({"reference", "query", ""});

	json configJson = CpuSwConfig();
	for (const auto& [gname, gvalues] : configJson.items()) {
		for (const auto& [optname, defaultValue] : gvalues.items()) {
			switch (defaultValue.type()) {
				case json::value_t::number_integer:
					options.add_option(gname, "", optname, "", cxxopts::value<int>(), std::to_string(defaultValue.get<int>()));
					break;
				case json::value_t::string:
					options.add_option(gname, "", optname, "", cxxopts::value<std::string>(), defaultValue);
					break;
				default:
					throw std::runtime_error("unsupported type");
					break;
			}
		}
	}

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

	for (const auto& [gname, gvalues] : configJson.items()) {
		for (const auto& [optName, optValue] : gvalues.items()) {
			if (result.count(optName)) {
				switch (optValue.type()) {
					case json::value_t::number_integer:
						configJson[gname][optName] = result[optName].as<int>();
						break;
					case json::value_t::string:
						configJson[gname][optName] = result[optName].as<std::string>();
						break;
					default:
						throw std::runtime_error("unsupported type");
						break;
				}
			}
		}
	}

	CpuSwConfig config = configJson.get<CpuSwConfig>();

	std::string refPath = result["reference"].as<std::string>();
	std::string queryPath = result["query"].as<std::string>();

	run_comparison(config, refPath, queryPath);
}