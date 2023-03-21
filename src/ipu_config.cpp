#include "ipu_config.h"

using json = nlohmann::json;

namespace ipu {
int get32bSize(int size8b) {
        return std::ceil(static_cast<double>(size8b) / 4.0);
}
int IPUAlgoConfig::getBufsize32b() const { return get32bSize(vertexBufferSize); }

int IPUAlgoConfig::getTotalNumberOfComparisons() const { return numVertices * maxComparisonsPerVertex; }

int IPUAlgoConfig::getMetaBufferSize32b() const { 
        return getTotalNumberOfComparisons() * getMetaStructSize32b();
};

bool IPUAlgoConfig::hasSeeds() const {
        return isSeeded(vtype);
}

int IPUAlgoConfig::getMetaStructSize32b() const {
        switch(this->vtype) {
        case VertexType::xdroprestrictedseedextend:
                return get32bSize(sizeof(XDropMeta));
                break;
        default:
                return get32bSize(sizeof(SWMeta));
                break;
        }
}

int IPUAlgoConfig::getVertexBufsize32b() const { return numVertices * getBufsize32b(); }

int IPUAlgoConfig::getInputBufferSize32b() const { return getVertexBufsize32b() + getMetaBufferSize32b(); }

size_t IPUAlgoConfig::getOffsetInputSequence() const { return 0; }

size_t IPUAlgoConfig::getOffsetMetadata() const {
        return getVertexBufsize32b();
}

void to_json(json& j, const SWConfig& c) {
        j = json{
                {"gapInit", c.gapInit},
                {"gapExtend", c.gapExtend},
                {"matchValue", c.matchValue},
                {"mismatchValue", c.mismatchValue},
                {"ambiguityValue", c.ambiguityValue},
                {"similarity", c.similarity},
                {"datatype", c.datatype},
                {"xDrop", c.xDrop},
                {"seedLength", c.seedLength}
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
        j.at("seedLength").get_to(c.seedLength);
        j.at("xDrop").get_to(c.xDrop);
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

void from_json(const json& j, Complexity& t) {
        t = ipu::strToComplexity(j);
}

void to_json(json& j, const Complexity& t) {
        j = ipu::complexityToConfigString(t);
}
}