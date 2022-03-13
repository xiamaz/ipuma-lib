#ifndef IPU_WORKER_HPP
#define IPU_WORKER_HPP
#include <plog/Log.h>
#include "ipumultidriver.hpp"
#include "batch.hpp"

void runIpuWorker(const int workerId, const int numWorkers, IPUMultiDriver& driver, const ipu::RawSequences& A, const ipu::RawSequences& B, const Batches& batches, swatlib::TickTock& inner) {
	std::vector<int32_t> inputBuffer(driver.config.ipuconfig.getInputBufferSize32b());
	std::vector<int32_t> resultBuffer(driver.config.ipuconfig.getTotalNumberOfComparisons() * 3);
	std::vector<int32_t> mappingBuffer(driver.config.ipuconfig.getTotalNumberOfComparisons(), 0);
	ipu::batchaffine::BlockAlignmentResults results;

	auto submit = [&](const Batch& batch) {
		std::vector<std::string> aBatch(batch.numCmps);
		std::vector<std::string> bBatch(batch.numCmps);

		for (int i = 0; i < batch.numCmps; ++i) {
			aBatch[i] = A.at(batch.startIndex + i);
			bBatch[i] = B.at(batch.startIndex + i);
		}
		driver.fill_input_buffer(aBatch, bBatch, inputBuffer, mappingBuffer);
		inner.tick();
		driver.submitWait(inputBuffer, resultBuffer);
		inner.tock();
		driver.fillResults(resultBuffer, mappingBuffer, results, batch.numCmps);
		driver.totalBatchesProcessed += 1;
		driver.totalCmpsProcessed += batch.numCmps;
	};

	int batchesProcessed = 0;
	for (int i = workerId; i < batches.size(); i += numWorkers) {
		submit(batches[i]);
		batchesProcessed++;
	}
	PLOGD.printf("Worker %d processed a total of %d batches.", workerId, batchesProcessed);
}

#endif