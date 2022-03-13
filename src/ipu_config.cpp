#include "ipu_config.h"

using json = nlohmann::json;

namespace ipu {
int IPUAlgoConfig::getBufsize32b() const { return std::ceil(static_cast<double>(bufsize) / 4.0); }

int IPUAlgoConfig::getTotalNumberOfComparisons() const { return tilesUsed * maxBatches; }

int IPUAlgoConfig::getMetaBufferSize32b() const { return getTotalNumberOfComparisons() * 4; };

int IPUAlgoConfig::getTotalBufsize32b() const { return tilesUsed * getBufsize32b(); }

int IPUAlgoConfig::getInputBufferSize32b() const { return getTotalBufsize32b() + getMetaBufferSize32b(); }

void to_json(json& j, const SWConfig& c) {
        j = json{
                {"gapInit", c.gapInit},
                {"gapExtend", c.gapExtend},
                {"matchValue", c.matchValue},
                {"mismatchValue", c.ambiguityValue},
                {"ambiguityValue", c.ambiguityValue},
                {"similarity", c.similarity},
                {"datatype", c.datatype},
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
                {"fillAlgo", c.fillAlgo},
                {"forwardOnly", c.forwardOnly},
                {"useRemoteBuffer", c.useRemoteBuffer},
                {"transmissionPrograms", c.transmissionPrograms},
                {"ioTiles", c.ioTiles}
        };
}

void from_json(const json& j, IPUAlgoConfig& c) {
        j.at("tilesUsed").get_to(c.tilesUsed);
        j.at("maxAB").get_to(c.maxAB);
        j.at("maxBatches").get_to(c.maxBatches);
        j.at("bufsize").get_to(c.bufsize);
        j.at("vtype").get_to(c.vtype);
        j.at("fillAlgo").get_to(c.fillAlgo);
        j.at("forwardOnly").get_to(c.forwardOnly);
        j.at("useRemoteBuffer").get_to(c.useRemoteBuffer);
        j.at("transmissionPrograms").get_to(c.transmissionPrograms);
        j.at("ioTiles").get_to(c.ioTiles);
}
}