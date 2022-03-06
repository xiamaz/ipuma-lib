#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "ipuswconfig.hpp"
#include "run_comparison.hpp"

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
		("reference", "Reference File (either fa or txt)", cxxopts::value<std::string>())
		("query", "Query File (either fa or txt)", cxxopts::value<std::string>())
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;

	options.positional_help("[reference_file] [query_file]");
	options.parse_positional({"reference", "query", ""});

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

	std::string refPath = result["reference"].as<std::string>();
	std::string queryPath = result["query"].as<std::string>();

	run_comparison(config, refPath, queryPath);
}