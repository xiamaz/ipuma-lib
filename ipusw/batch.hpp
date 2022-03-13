#ifndef BATCH_HPP
#define BATCH_HPP
#include <stdexcept>
#include <vector>
#include <string>

#include <plog/Log.h>

struct Batch {
	int startIndex;
	int endIndex;
	int numCmps;
	int totalSize;
};

using Batches = std::vector<Batch>;

Batches createBatches(std::vector<std::string>& A, std::vector<std::string>& B, const int batchCmpLimit, const int batchDataLimit) {
	if (A.size() != B.size()) {
		throw std::runtime_error("Sizes of A and B are not equal");
	}
	Batches batches;
	int numCmps = 0;
	int totalSize = 0;
	int batchBegin = 0;
	for (int i = 0; i < A.size(); ++i) {
		const auto alen = A[i].size();
		const auto blen = B[i].size();

		if ((numCmps + 1 < batchCmpLimit) && (totalSize + alen + blen < batchDataLimit)) {
			numCmps++;
			totalSize += alen + blen;
		} else {
			batches.push_back({.startIndex = batchBegin, .endIndex = batchBegin + numCmps, .numCmps = numCmps, .totalSize = totalSize});

			batchBegin += numCmps;
			numCmps = 1;
			totalSize = alen + blen;
		}
	}
	if (numCmps > 0) {
		batches.push_back({.startIndex = batchBegin, .endIndex = batchBegin + numCmps, .numCmps = numCmps, .totalSize = totalSize});
	}

	PLOGD.printf("Created %d batches for a total of %d comparisons.", batches.size(), A.size());
	return batches;
}

#endif