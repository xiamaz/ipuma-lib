#include "types.h"

#include <sstream>
#include <stdexcept>
#include <string>

namespace ipu {
std::string vertexTypeToIpuLabel(VertexType v) { return typeLabels[static_cast<int>(v)]; }
std::string vertexTypeToConfigString(VertexType v) { return vertexTypeNames[static_cast<int>(v)]; }
std::string algorithmToConfigString(Algorithm a) { return algoNames[static_cast<int>(a)]; }
std::string complexityToConfigString(Complexity c) { return complexityNames[static_cast<int>(c)]; }

VertexType strToVertexType(std::string s) {
  for (int i = 0; i < vertexTypeNames.size(); ++i) {
    if (vertexTypeNames[i] == s) return static_cast<VertexType>(i);
  }
  throw std::runtime_error("Invalid vertex type string: " + s);
}

Algorithm strToAlgorithm(std::string s) {
  for (int i = 0; i < algoNames.size(); ++i) {
    if (algoNames[i] == s) return static_cast<Algorithm>(i);
  }
  throw std::runtime_error("Invalid algorithm string: " + s);
}

Complexity strToComplexity(std::string s) {
  for (int i = 0; i < complexityNames.size(); ++i) {
    if (complexityNames[i] == s) return static_cast<Complexity>(i);
  }
  throw std::runtime_error("Invalid complexity string: " + s);
}

template<>
void add_comparison(std::vector<Comparison>& cmps, Comparison& cmp, int seedLen) {
	cmps.push_back(cmp);
}

template<>
void add_comparison(std::vector<MultiComparison>& cmps, Comparison& cmp, int seedLen) {
	cmps.push_back({{cmp}, seedLen});
}

bool isSeeded(VertexType vtype) {
  return vtype == VertexType::xdroprestrictedseedextend;
}

bool Comparison::operator<(const Comparison& other) const {
  return this->complexity < other.complexity;
}

bool MultiComparison::operator<(const MultiComparison& other) const {
  return this->complexity < other.complexity;
}

Comparisons convertToComparisons(const MultiComparisons& mcmps) {
  Comparisons cmps;

  for (const auto& m : mcmps) {
    for (const auto& c : m.comparisons) {
      cmps.push_back(c);
    }
  }
  return std::move(cmps);
}

std::string Comparison::toString() const {
  std::stringstream ss;
  ss << "Comp[" << originalComparisonIndex << ": a(l" << sizeA << ") b(l" << sizeB << ")]";
  return ss.str();
}

MultiComparison::MultiComparison() {}

MultiComparison::MultiComparison(const std::vector<Comparison>& cmps, const int seedLength) : comparisons(cmps), comparisonCount(cmps.size() * NSEEDS) {
  for (const auto& comparison : comparisons) {
    for (const auto& pair : comparison.seeds) {
      int left = pair.seedAStartPos * pair.seedBStartPos;
      int right = (comparison.sizeA - seedLength - pair.seedAStartPos) * (comparison.sizeB - seedLength - pair.seedBStartPos);
      complexity += (pair.seedAStartPos != -1 ? 1 : 0) * (left + right);
    }
    seqs[comparison.indexA] = comparison.sizeA;
    seqs[comparison.indexB] = comparison.sizeB;
  }
  for (const auto& [key, value] : seqs) {
    totalSeqSize += value;
  }
}

void to_json(json& j, const Comparison& c) {
  j = json{
      {"originalComparisonIndex", c.originalComparisonIndex},
      {"indexA", c.indexA},
      {"sizeA", c.sizeA},
      {"indexB", c.indexB},
      {"sizeB", c.sizeB},
      {"complexity", c.complexity},
      {"seeds", c.seeds},
  };
}

void from_json(const json& j, Comparison& c) {
	j.at("originalComparisonIndex").get_to(c.originalComparisonIndex);
	j.at("indexA").get_to(c.indexA);
	j.at("indexB").get_to(c.indexB);
	j.at("sizeA").get_to(c.sizeA);
	j.at("sizeB").get_to(c.sizeB);
	j.at("complexity").get_to(c.complexity);
	j.at("seeds").get_to(c.seeds);
}

void to_json(json& j, const MultiComparison& c) {
  j = json{
      {"totalSeqSize", c.totalSeqSize},
      {"complexity", c.complexity},
      {"comparisonCount", c.comparisonCount},
      {"comparisons", c.comparisons},
      {"seqs", c.seqs},
  };
}

void from_json(const json& j, MultiComparison& c) {
	j.at("totalSeqSize").get_to(c.totalSeqSize);
	j.at("complexity").get_to(c.complexity);
	j.at("comparisonCount").get_to(c.comparisonCount);
	j.at("comparisons").get_to(c.comparisons);
	j.at("seqs").get_to(c.seqs);
}

void to_json(json& j, const SeedPair& c) {
	j = { c.seedAStartPos, c.seedBStartPos };
}
void from_json(const json& j, SeedPair& c) {
	int32_t aPos;
	int32_t bPos;
	j.at(0).get_to(aPos);
	j.at(1).get_to(bPos);
	c.seedAStartPos = aPos;
	c.seedBStartPos = bPos;
}
}  // namespace ipu