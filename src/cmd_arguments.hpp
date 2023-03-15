#pragma once
#include <nlohmann/json.hpp>
#include <cxxopts.hpp>

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
			throw std::runtime_error("unsupported type: " + std::string(v.type_name()) + " for field: " + k);
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
				throw std::runtime_error("unsupported type: " + std::string(v.type_name()) + " for field: " + k);
				break;
		}
	}
}