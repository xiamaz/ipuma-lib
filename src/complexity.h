#pragma once

#include "types.h"

namespace ipu {
	size_t calculateComplexity(const Comparison& cmp, Complexity algo);
	size_t calculateComplexity(const MultiComparison& cmp, Complexity algo);

	void addComplexity(Comparisons& cmps, Complexity algo);
	void addComplexity(MultiComparisons& cmps, Complexity algo);
}