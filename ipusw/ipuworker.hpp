#ifndef IPU_WORKER_HPP
#define IPU_WORKER_HPP
#include <plog/Log.h>
#include "ipumultidriver.hpp"
#include "batch.hpp"

// worker for getting results
void workerResult(const int workerId, const int numWorkers, IPUMultiDriver& driver, Batches& batches) {
	int finished = 0;
	int workerBatchCount = 0;
	for (int i = workerId; i < batches.size(); i += numWorkers) {
		workerBatchCount++;
	}
	while (finished < workerBatchCount) {
		for (int i = workerId; i < batches.size(); i += numWorkers) {
			auto& batch = batches[i];
			if (batch.received) continue;
			if (batch.job != nullptr) {
				driver.wait(batch.job);
				ipu::batchaffine::BlockAlignmentResults results;
		    driver.fillResults(batch.resultBuffer, batch.mappingBuffer, results, batch.numCmps);
				batch.received = true;
				finished++;
			}
		}
	}
	PLOGW << "Joining result worker";
}

void runIpuWorker(const int workerId, const int numWorkers, IPUMultiDriver& driver, const ipu::RawSequences& A, const ipu::RawSequences& B, Batches& batches) {
	const auto inputBufferSize = driver.config.ipuconfig.getInputBufferSize32b();
	const auto mappingBufferSize = driver.config.ipuconfig.getTotalNumberOfComparisons();
	const auto resultBufferSize = driver.config.ipuconfig.getTotalNumberOfComparisons() * 3;
	// std::vector<int32_t> inputBuffer(driver.config.ipuconfig.getInputBufferSize32b());
	// std::vector<int32_t> resultBuffer(driver.config.ipuconfig.getTotalNumberOfComparisons() * 3);
	// std::vector<int32_t> mappingBuffer(driver.config.ipuconfig.getTotalNumberOfComparisons(), 0);
	// ipu::batchaffine::BlockAlignmentResults results;

	auto submit = [&](Batch& batch) {
		std::vector<std::string> aBatch(batch.numCmps);
		std::vector<std::string> bBatch(batch.numCmps);

		for (int i = 0; i < batch.numCmps; ++i) {
			aBatch[i] = A.at(batch.startIndex + i);
			bBatch[i] = B.at(batch.startIndex + i);
		}
		batch.inputBuffer.resize(inputBufferSize);
		batch.mappingBuffer.resize(mappingBufferSize);
		batch.resultBuffer.resize(resultBufferSize);
		driver.fill_input_buffer(aBatch, bBatch, batch.inputBuffer, batch.mappingBuffer);
		// inner.tick();
		// driver.submitWait(inputBuffer, resultBuffer);
		// inner.tock();
		// driver.fillResults(resultBuffer, mappingBuffer, results, batch.numCmps);
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