#ifndef IPU_SW_CONFIG_HPP
#define IPU_SW_CONFIG_HPP

#include "swatlib/swatlib.h"
#include "driver.hpp"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace swatlib {
	void to_json(json& j, const swatlib::DataType& d) {
		j = swatlib::dataTypeToStr(d);
	}

	void from_json (const json& j, swatlib::DataType& d) {
		d = swatlib::strToDataType(j);
	}

	void to_json(json& j, const swatlib::Similarity& s) {
		j = swatlib::similarityToStr(s);
	}

	void from_json(const json& j, swatlib::Similarity& s) {
		s = swatlib::strToSimilarity(j);
	}
}

namespace ipu {
	void to_json(json& j, const SWConfig& c) {
		j = json{
			{"gapInit", c.gapInit},
			{"gapExtend", c.gapExtend},
			{"matchValue", c.matchValue},
			{"mismatchValue", c.ambiguityValue},
			{"ambiguityValue", c.ambiguityValue},
			{"similarity", c.similarity},
			{"datatype", c.datatype}
		};
	}

	void from_json(const json& j, SWConfig& c) {
		j.at("gapInit").get_to(c.gapInit);
		j.at("gapExtend").get_to(c.gapExtend);
		j.at("matchValue").get_to(c.matchValue);
		j.at("mismatchValue").get_to(c.mismatchValue);
		j.at("ambiguityValue").get_to(c.ambiguityValue);
		j.at("similarity").get_to(c.similarity);
		j.at("datatype").get_to(c.datatype);
	}

	void from_json(const json& j, Algorithm& a) {
		a = ipu::strToAlgorithm(j);
	}

	void to_json(json& j, const Algorithm& a) {
		j = ipu::algorithmToConfigString(a);
	}

	void from_json(const json& j, VertexType& t) {
		t = ipu::strToVertexType(j);
	}

	void to_json(json& j, const VertexType& t) {
		j = ipu::vertexTypeToConfigString(t);
	}

	void to_json(json& j, const IPUAlgoConfig& c) {
		j = json{
			{"tilesUsed", c.tilesUsed},
			{"maxAB", c.maxAB},
			{"maxBatches", c.maxBatches},
			{"bufsize", c.bufsize},
			{"vtype", c.vtype},
			{"fillAlgo", c.fillAlgo}
		};
	}

	void from_json(const json& j, IPUAlgoConfig& c) {
		j.at("tilesUsed").get_to(c.tilesUsed);
		j.at("maxAB").get_to(c.maxAB);
		j.at("maxBatches").get_to(c.maxBatches);
		j.at("bufsize").get_to(c.bufsize);
		j.at("vtype").get_to(c.vtype);
		j.at("fillAlgo").get_to(c.fillAlgo);
	}
}

struct IpuSwConfig {
	ipu::SWConfig swconfig;
	ipu::IPUAlgoConfig ipuconfig;

	IpuSwConfig() : swconfig(SW_CONFIGURATION), ipuconfig(ALGO_CONFIGURATION) {}

	IpuSwConfig(ipu::SWConfig sw, ipu::IPUAlgoConfig ipu) : swconfig(sw), ipuconfig(ipu) {}
};

void to_json(json& j, const IpuSwConfig& c) {
	j = json{
		{"sw", c.swconfig},
		{"ipu", c.ipuconfig}
	};
}

void from_json(const json& j, IpuSwConfig& c) {
	j.at("sw").get_to(c.swconfig);
	j.at("ipu").get_to(c.ipuconfig);
}
#endif