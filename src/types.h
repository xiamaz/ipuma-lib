#ifndef IPUMA_TYPES_H
#define IPUMA_TYPES_H
#include <string>
#include <vector>

#define IPU_JSON_LOG_TAG "IPUSWLOG"

namespace ipu {
typedef std::vector<std::string> RawSequences;
typedef std::vector<uint8_t> EncSequences;

struct __attribute__((__packed__)) Comparison {
	int32_t indexA;
	int32_t indexB;
	int32_t seedAStartPos;	
	int32_t seedBStartPos;	
};

struct __attribute__((__packed__)) SWMeta {
	int32_t sizeA;
	int32_t offsetA;
	int32_t sizeB;
	int32_t offsetB;
};

struct __attribute__((__packed__)) XDropMeta : public SWMeta {
	int32_t seedAStartPos;	
	int32_t seedBStartPos;	
};

typedef std::vector<Comparison> Comparisons;

enum class VertexType { cpp, assembly, multi, multiasm, xdrop, multixdrop, greedyxdrop, multibandxdrop, xdropseedextend};
static const std::vector<std::string> vertexTypeNames = {"cpp", "assembly", "multi", "multiasm", "xdrop", "multixdrop", "greedyxdrop", "multibandxdrop", "xdropseedextend"};
static const std::string typeLabels[] = {"SWAffine", "SWAffineAsm", "MultiSWAffine", "MultiSWAffineAsm", "XDrop", "MultiXDrop", "GreedyXDrop", "MultiBandXDrop", "SeedExtendXDrop"};

enum class Algorithm {fillFirst, roundRobin, greedy};
static const std::vector<std::string> algoNames = {"fillfirst", "roundrobin", "greedy"};
std::string vertexTypeToIpuLabel(VertexType v);
std::string vertexTypeToConfigString(VertexType v);
std::string algorithmToConfigString(Algorithm a);
VertexType strToVertexType(std::string s);
Algorithm strToAlgorithm(std::string s);
}

#endif