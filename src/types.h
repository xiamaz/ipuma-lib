#ifndef IPUMA_TYPES_H
#define IPUMA_TYPES_H
#include <string>
#include <vector>

namespace ipu {
typedef std::vector<std::string> RawSequences;
typedef std::vector<uint8_t> EncSequences;

struct __attribute__((__packed__)) Comparison {
	int32_t indexA;
	int32_t indexB;
};

typedef std::vector<Comparison> Comparisons;

enum class VertexType { cpp, assembly, multi, multiasm, stripedasm, multistriped, multistripedasm };
static const std::vector<std::string> vertexTypeNames = {"cpp", "assembly", "multi", "multiasm", "stripedasm", "multistriped", "multistripedasm"};
enum class Algorithm {fillFirst, roundRobin, greedy};
static const std::vector<std::string> algoNames = {"fillfirst", "roundrobin", "greedy"};

static const std::string typeLabels[] = {"SWAffine", "SWAffineAsm", "MultiSWAffine", "MultiSWAffineAsm", "StripedSWAffineAsm", "MultiSWStriped", "MultiSWStripedAsm"};
std::string vertexTypeToIpuLabel(VertexType v);
std::string vertexTypeToConfigString(VertexType v);
std::string algorithmToConfigString(Algorithm a);
VertexType strToVertexType(std::string s);
Algorithm strToAlgorithm(std::string s);
}

#endif