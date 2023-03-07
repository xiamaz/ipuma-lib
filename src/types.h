#ifndef IPUMA_TYPES_H
#define IPUMA_TYPES_H
#include <string>
#include <string_view>
#include <vector>
#include <array>
#include "types.h"
#include "shared_types.h"

#define IPU_JSON_LOG_TAG "IPUSWLOG"

namespace ipu {
struct Comparison {
	int32_t indexA;
	int32_t indexB;
	std::array<SeedPair, NSEEDS> seeds;
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