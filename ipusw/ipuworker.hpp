#ifndef IPU_WORKER_HPP
#define IPU_WORKER_HPP
#include <plog/Log.h>
#include "ipumultidriver.hpp"

void runIpuWorker(const int workerId, const int numWorkers, IPUMultiDriver& driver, const ipu::RawSequences& A, const ipu::RawSequences& B, swatlib::TickTock& inner) {
	std::vector<int32_t> inputBuffer(driver.config.ipuconfig.getInputBufferSize32b());
	std::vector<int32_t> resultBuffer(driver.config.ipuconfig.getTotalNumberOfComparisons() * 3);
	std::vector<int32_t> mappingBuffer(driver.config.ipuconfig.getTotalNumberOfComparisons(), 0);
	ipu::batchaffine::BlockAlignmentResults results;

  int seqsPerWorker = (A.size() + numWorkers - 1) / numWorkers;
  int startIndex = std::min(seqsPerWorker * workerId, static_cast<int>(A.size()));
  int endIndex = std::min(seqsPerWorker * (workerId + 1), static_cast<int>(A.size()));
	PLOGD.printf("Worker %d: processing %d sequences from %d to %d (actual %d)", workerId, seqsPerWorker, startIndex, endIndex, endIndex - startIndex);

	auto submit = [&](int batchBegin, int numCmps) {
		std::vector<std::string> aBatch(numCmps);
		std::vector<std::string> bBatch(numCmps);

		for (int i = 0; i < numCmps; ++i)
		{
			aBatch[i] = A.at(batchBegin + i);
			bBatch[i] = B.at(batchBegin + i);
		}
		driver.fill_input_buffer(aBatch, bBatch, inputBuffer, mappingBuffer);
		inner.tick();
		driver.submitWait(inputBuffer, resultBuffer);
		inner.tock();
		driver.fillResults(resultBuffer, mappingBuffer, results, numCmps);
	};

  int batchCmpLimit = driver.config.ipuconfig.getTotalNumberOfComparisons() - driver.config.ipuconfig.maxBatches;
  int slack = driver.config.ipuconfig.bufsize * 100; // only fill buffer up to 90%, hack to make sure we don't fail
  int batchDataLimit = (driver.config.ipuconfig.getTotalBufsize32b() * 4) - slack;

	int numCmps = 0;
	int totalSize = 0;
	int batchBegin = startIndex;

	for (int i = startIndex; i < endIndex; ++i)
	{
		const auto alen = A[i].size();
    const auto blen = B[i].size();
		if ((numCmps + 1 < batchCmpLimit) && (totalSize + alen + blen < batchDataLimit))
		{
			numCmps++;
			totalSize += alen + blen;
		}
		else
		{
			submit(batchBegin, numCmps);

			batchBegin += numCmps;
			numCmps = 1;
			totalSize = alen + blen;
		}
	}
	if (numCmps > 0)
	{
		submit(batchBegin, numCmps);
	}

}

#endif