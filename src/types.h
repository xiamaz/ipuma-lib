#ifndef IPUMA_TYPES_H
#define IPUMA_TYPES_H
#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <unordered_map>
#include "types.h"
#include "shared_types.h"

#include "nlohmann/json.hpp"

#define IPU_JSON_LOG_TAG "IPUSWLOG"

namespace ipu {

using json = nlohmann::json;

struct Comparison {
	int64_t originalComparisonIndex;
	int32_t indexA;
	int32_t sizeA;
	int32_t indexB;
	int32_t sizeB;
	std::array<SeedPair, NSEEDS> seeds;
	size_t complexity;

  bool operator<(const Comparison& other) const;

	std::string toString() const;
};

struct MultiComparison {
	int32_t totalSeqSize = 0;
	int32_t complexity = 0;
	int32_t comparisonCount = 0;
	std::vector<Comparison> comparisons;
	std::unordered_map<int32_t, int32_t> seqs; // index -> sequence_length

  bool operator<(const MultiComparison& other) const;

	MultiComparison();

	MultiComparison(const std::vector<Comparison>& cmps, const int seedLength);
};

typedef std::vector<std::string_view> RawSequences;
typedef std::vector<uint8_t> EncSequences;


typedef std::vector<Comparison> Comparisons;
typedef std::vector<MultiComparison> MultiComparisons;

template<typename C>
void add_comparison(std::vector<C>& cmps, Comparison& cmp, int seedLen);

Comparisons convertToComparisons(const MultiComparisons&);

enum class VertexType { cpp, assembly, multi, multiasm, xdroprestrictedseedextend};
static const std::vector<std::string> vertexTypeNames = {"cpp", "assembly", "multi", "multiasm", "xdroprestrictedseedextend"};
static const std::string typeLabels[] = {"SWAffine", "SWAffineAsm", "MultiSWAffine", "MultiSWAffineAsm", "SeedExtendRestrictedXDrop"};

bool isSeeded(VertexType vtype);

enum class Algorithm {fillFirst, roundRobin, greedy};
static const std::vector<std::string> algoNames = {"fillfirst", "roundrobin", "greedy"};

enum class Complexity {precomputed, cellcount, sequence_length, xdrop};
static const std::vector<std::string> complexityNames = {"precomputed", "cellcount", "sequence_length", "xdrop"};

enum class PartitionAdd {sequential, alternating, heap};
NLOHMANN_JSON_SERIALIZE_ENUM(PartitionAdd, {
    {PartitionAdd::sequential, "sequential"},
    {PartitionAdd::alternating, "alternating"},
    {PartitionAdd::heap, "heap"},
});


std::string vertexTypeToIpuLabel(VertexType v);
std::string vertexTypeToConfigString(VertexType v);
VertexType strToVertexType(std::string s);

std::string algorithmToConfigString(Algorithm a);
Algorithm strToAlgorithm(std::string s);

std::string complexityToConfigString(Complexity);
Complexity strToComplexity(std::string s);

void to_json(json& j, const Comparison& c);
void from_json(const json& j, Comparison& c);
void to_json(json& j, const MultiComparison& c);
void from_json(const json& j, MultiComparison& c);
void to_json(json& j, const SeedPair& c);
void from_json(const json& j, SeedPair& c);
}

#endif