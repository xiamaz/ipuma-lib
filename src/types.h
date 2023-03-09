#ifndef IPUMA_TYPES_H
#define IPUMA_TYPES_H
#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <unordered_map>
#include "types.h"
#include "shared_types.h"

#define IPU_JSON_LOG_TAG "IPUSWLOG"

namespace ipu {

struct Comparison {
	int64_t originalComparisonIndex;
	int32_t indexA;
	int32_t sizeA;
	int32_t indexB;
	int32_t sizeB;
	std::array<SeedPair, NSEEDS> seeds;
	// int real_complexity = 0;
};

struct MultiComparison {
private:
	int32_t totalSeqSize = 0;
	int32_t complexity = 0;
	int32_t comparisonCount = 0;
	std::unordered_map<int32_t, int32_t> seqs;
public:
	std::vector<Comparison> comparisons;
	MultiComparison(const std::vector<Comparison>& cmps, const int seedLength) : comparisons(cmps), comparisonCount(cmps.size() * NSEEDS) {
		for (const auto& comparison : comparisons) {
    			for (const auto &pair : comparison.seeds) {
    				int left = pair.seedAStartPos * pair.seedBStartPos;
    				int right = (comparison.sizeA - seedLength - pair.seedAStartPos) * (comparison.sizeB - seedLength - pair.seedBStartPos);
    				complexity += (pair.seedAStartPos != -1 ? 1 : 0) * (left + right);
    			}
			seqs[comparison.indexA] = comparison.sizeA;
			seqs[comparison.indexB] = comparison.sizeB;
		}
		for (const auto &[key, value]: seqs) {
			totalSeqSize += value;
		}
	}
};

typedef std::vector<std::string_view> RawSequences;
typedef std::vector<uint8_t> EncSequences;


typedef std::vector<Comparison> Comparisons;

enum class VertexType { cpp, assembly, multi, multiasm, xdrop, multixdrop, greedyxdrop, multibandxdrop, xdropseedextend, xdroprestrictedseedextend};
static const std::vector<std::string> vertexTypeNames = {"cpp", "assembly", "multi", "multiasm", "xdrop", "multixdrop", "greedyxdrop", "multibandxdrop", "xdropseedextend", "xdroprestrictedseedextend"};
static const std::string typeLabels[] = {"SWAffine", "SWAffineAsm", "MultiSWAffine", "MultiSWAffineAsm", "XDrop", "MultiXDrop", "GreedyXDrop", "MultiBandXDrop", "SeedExtendXDrop", "SeedExtendRestrictedXDrop"};

enum class Algorithm {fillFirst, roundRobin, greedy};
static const std::vector<std::string> algoNames = {"fillfirst", "roundrobin", "greedy"};
std::string vertexTypeToIpuLabel(VertexType v);
std::string vertexTypeToConfigString(VertexType v);
std::string algorithmToConfigString(Algorithm a);
VertexType strToVertexType(std::string s);
Algorithm strToAlgorithm(std::string s);
}

#endif