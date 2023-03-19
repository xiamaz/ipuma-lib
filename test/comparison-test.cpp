#include "gtest/gtest.h"
#include <plog/Log.h>

#include "types.h"

using namespace ipu;
using json = nlohmann::json;

void compareComparison(Comparison& c, Comparison& cj) {
	ASSERT_EQ(cj.originalComparisonIndex, c.originalComparisonIndex);
	ASSERT_EQ(cj.indexA, c.indexA);
	ASSERT_EQ(cj.sizeA, c.sizeA);
	ASSERT_EQ(cj.indexB, c.indexB);
	ASSERT_EQ(cj.sizeB, c.sizeB);
	ASSERT_EQ(cj.seeds.size(), c.seeds.size());
	ASSERT_EQ(cj.complexity, c.complexity);
	for (int i = 0; i < cj.seeds.size(); ++i) {
		ASSERT_EQ(cj.seeds[i].seedAStartPos, c.seeds[i].seedAStartPos);
		ASSERT_EQ(cj.seeds[i].seedBStartPos, c.seeds[i].seedBStartPos);
	}
}

TEST(ComparisonTest, testJsonSerialize) {
	Comparison c{
		.originalComparisonIndex = 1,
		.indexA = 2,
		.sizeA = 3,
		.indexB = 4,
		.sizeB = 5,
		.seeds = {{{10, 12}}},
		.complexity = 42
	};

	json j = c;
	Comparison cj = j;

	compareComparison(c, cj);
}

TEST(MultiComparisonTest, testJsonSerialize) {
	Comparison c{
		.originalComparisonIndex = 1,
		.indexA = 2,
		.sizeA = 3,
		.indexB = 4,
		.sizeB = 5,
		.seeds = {{{10, 12}}},
		.complexity = 42
	};

	MultiComparison mc({c}, 17);

	json j = mc;
	MultiComparison mcj = j;
	ASSERT_EQ(mcj.totalSeqSize, mc.totalSeqSize);
	ASSERT_EQ(mcj.complexity, mc.complexity);
	ASSERT_EQ(mcj.comparisonCount, mc.comparisonCount);
	ASSERT_EQ(mcj.comparisons.size(), mc.comparisons.size());
	compareComparison(mcj.comparisons[0], mc.comparisons[0]);
}