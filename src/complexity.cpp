#include "complexity.h"

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
	}
}