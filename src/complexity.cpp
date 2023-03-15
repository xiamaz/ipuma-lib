#include "complexity.h"
#include <stdexcept>
#include <plog/Log.h>

namespace ipu {
	size_t xDropComplexity(const Comparison& cmp) {
		size_t complexity = 0;
    for (const auto &pair : cmp.seeds) {
     int left = pair.seedAStartPos * pair.seedBStartPos;
     int right = (cmp.sizeA - 17 - pair.seedAStartPos) * (cmp.sizeB - 17 - pair.seedBStartPos);
     complexity += (pair.seedAStartPos != -1 ? 1 : 0) * (left + right);
    }
		return complexity;
	}

	size_t calculateComplexity(const MultiComparison& cmp, Complexity algo) {
		size_t cellCount = 0;
		size_t totalComplexity = 0;
		switch (algo) {
		case Complexity::precomputed:
			return cmp.complexity;
			break;
		case Complexity::sequence_length:
			return cmp.totalSeqSize;
			break;
		case Complexity::cellcount:
			for (const auto& c : cmp.comparisons) {
				cellCount += c.sizeA * c.sizeB * NSEEDS;
			}
			return cellCount;
			break;
		case Complexity::xdrop:
			for (const auto& c : cmp.comparisons) {
				totalComplexity += xDropComplexity(c);
			}
			return totalComplexity;
			break;
		}
		throw std::runtime_error("Unhandled complexity algorithm type.");
	}

	size_t calculateComplexity(const Comparison& cmp, Complexity algo) {
		switch (algo) {
		case Complexity::precomputed:
			return cmp.complexity;
			break;
		case Complexity::sequence_length:
			return cmp.sizeA + cmp.sizeB;
			break;
		case Complexity::cellcount:
			return cmp.sizeA * cmp.sizeB * NSEEDS;
			break;
		case Complexity::xdrop:
			return xDropComplexity(cmp);
			break;
		}
		throw std::runtime_error("Unhandled complexity algorithm type.");
	}

	void addComplexity(Comparisons& Cmps, Complexity algo) {
    PLOGD << "Set complexity using " << complexityToConfigString(algo);
    for (auto&& cmp : Cmps) {
      cmp.complexity = calculateComplexity(cmp, algo);
    }
	}

	void addComplexity(MultiComparisons& Cmps, Complexity algo) {
    PLOGD << "Set complexity using " << complexityToConfigString(algo);
    for (auto&& cmp : Cmps) {
      cmp.complexity = calculateComplexity(cmp, algo);
    }
	}
}