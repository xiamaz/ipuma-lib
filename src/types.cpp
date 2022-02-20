#include "types.h"
#include <string>
#include <stdexcept>

namespace ipu {
std::string vertexTypeToIpuLabel(VertexType v) { return typeLabels[static_cast<int>(v)]; }
std::string vertexTypeToConfigString(VertexType v) { return vertexTypeNames[static_cast<int>(v)]; }
std::string algorithmToConfigString(Algorithm a) { return algoNames[static_cast<int>(a)]; }
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
}