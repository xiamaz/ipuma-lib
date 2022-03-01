#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include "ipuswconfig.hpp"
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <omp.h>

#include <plog/Log.h>

using json = nlohmann::json;

using Fastas = std::vector<swatlib::Fasta>;

void load_data(const std::string& path, std::vector<std::string>& seqs, int count = 0) {
	if (count > 0) {
		PLOGI << "Loading " << count << " entries from " << path;
	} else {
		PLOGI << "Loading all entries from " << path;
	}
  std::ifstream is;
  is.open(path);
  if (is.fail()) {
      throw std::runtime_error("Opening file at " + path + " failed.");
  }
  swatlib::Fasta f;
	int i = 0;
	while (!(count) || i < count) {
      if (is >> f) {
          seqs.push_back(f.sequence);
      } else {
          seqs.push_back(f.sequence);
          break;
      }
			++i;
	}
}

std::mutex ipuMutex;
std::mutex bufferSelectionMutex;
std::mutex prepareWorkerMutex;

int submittedBatches = 0;
int processedBatches = 0;
std::vector<int> bufferMaxBucket;

int startedWorkers = 0;
int stoppedWorkers = 0;
std::mutex runningWorkersMutex;

void ipu_init(const IpuSwConfig& config, std::vector<ipu::batchaffine::SWAlgorithm>& drivers) {
	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, true);
	ipuMutex.lock();
	drivers.push_back(std::move(driver));
	ipuMutex.unlock();
}

void ipu_prepare(const ipu::SWConfig& swconfig, const ipu::IPUAlgoConfig& ipuconfig, const int startIndex, const int endIndex, const ipu::RawSequences& A, const ipu::RawSequences& B, std::vector<int>& bufferCmps, std::vector<std::vector<int32_t>>& input_bufs, std::vector<std::vector<int>>& mapping_bufs) {
	runningWorkersMutex.lock();
	startedWorkers += 1;
	runningWorkersMutex.unlock();
	// heuritic slicing for greater speed
	int numCmps = 0;
	int totalSize = 0;

	ipu::RawSequences::const_iterator aBegin, aEnd, bBegin, bEnd;

	int batchCmpLimit = ipuconfig.getTotalNumberOfComparisons() - ipuconfig.maxBatches;
	int batchDataLimit = (ipuconfig.getTotalBufsize32b() * 4) - ipuconfig.bufsize;

	aBegin = A.begin() + startIndex;
	bBegin = B.begin() + startIndex;

	for (int i = startIndex; i < A.size() && i < endIndex; ++i) {
		const auto alen = A[i].size();
		const auto blen = B[i].size();
		if (numCmps + 1 < batchCmpLimit && totalSize + alen + blen < batchDataLimit) {
			numCmps++;
			totalSize += alen + blen;
		} else {
			aEnd = aBegin + numCmps;
			bEnd = bBegin + numCmps;

			std::vector<std::string> aBatch(aBegin, aEnd);
			std::vector<std::string> bBatch(bBegin, bEnd);

			int c;
			while (true) {
				bufferSelectionMutex.lock();
				for (c = 0; c < bufferCmps.size(); ++c) {
					if (bufferCmps[c] == 0) {
						bufferCmps[c] = -1; // reserve buffer
						break;
					}
				}
				bufferSelectionMutex.unlock();
				if (c < bufferCmps.size()) break;
			}
			auto maxBucket = ipu::batchaffine::SWAlgorithm::prepare_remote(swconfig, ipuconfig, aBatch, bBatch, &*input_bufs[c].begin(), &*input_bufs[c].end(), mapping_bufs[c].data());
			bufferSelectionMutex.lock();
			bufferCmps[c] = numCmps;
			bufferMaxBucket[c] = maxBucket;
			submittedBatches += 1;
			bufferSelectionMutex.unlock();
			PLOGI << "Putting " << numCmps << " into " << c;

			// trigger compare

			numCmps = 1;
			totalSize = alen + blen;
			aBegin = aEnd + 1;
			bBegin = bEnd + 1;
		}
	}
	if (numCmps > 0) {
		aEnd = aBegin + numCmps;
		bEnd = bBegin + numCmps;

		std::vector<std::string> aBatch(aBegin, aEnd);
		std::vector<std::string> bBatch(bBegin, bEnd);
		int c;
		while (true) {
			bufferSelectionMutex.lock();
			for (c = 0; c < bufferCmps.size(); ++c) {
				if (bufferCmps[c] == 0) {
					bufferCmps[c] = -1;
					break;
				}
			}
			bufferSelectionMutex.unlock();
			if (c < bufferCmps.size()) break;
		}
		auto maxBucket = ipu::batchaffine::SWAlgorithm::prepare_remote(swconfig, ipuconfig, aBatch, bBatch, &*input_bufs[c].begin(), &*input_bufs[c].end(), mapping_bufs[c].data());
		bufferSelectionMutex.lock();
		bufferCmps[c] = numCmps;
		bufferMaxBucket[c] = maxBucket;
		submittedBatches += 1;
		bufferSelectionMutex.unlock();
		PLOGI << "Putting " << numCmps << " into " << c;
		// trigger compare
	}

	PLOGW << "Shutting down";

	runningWorkersMutex.lock();
	stoppedWorkers += 1;
	runningWorkersMutex.unlock();
}

void ipu_run(ipu::batchaffine::SWAlgorithm& driver, const ipu::RawSequences& A, const ipu::RawSequences& B, const int workerId, const int numWorkers, swatlib::TickTock& t) {
	const auto& ipuconfig = driver.algoconfig;
	const auto& swconfig = driver.config;

	const int numInputBufs = 2;
	const int numPrepareWorkers = 1;

  std::vector<std::vector<int>> mapping_bufs(numInputBufs, std::vector<int>(ipuconfig.getTotalNumberOfComparisons(), 0));
  std::vector<std::vector<int32_t>> input_bufs(numInputBufs, std::vector<int32_t>(ipuconfig.getInputBufferSize32b()));
  size_t results_size = ipuconfig.getTotalNumberOfComparisons() * 3;
  std::vector<int32_t> results_buf(results_size);

	std::vector<int32_t> scores(A.size());
	std::vector<int32_t> arange(A.size());
	std::vector<int32_t> brange(A.size());

	int results_offset = 0;

	std::vector<int> bufferCmps(numInputBufs, 0);
	bufferMaxBucket.resize(numInputBufs);
	int seqsPerWorker = (A.size() + numWorkers - 1) / numWorkers;
	int startIndex = seqsPerWorker * workerId;
	int endIndex = seqsPerWorker * (workerId + 1);

	int workRange = endIndex - startIndex;
	int workChunks = (workRange + numPrepareWorkers - 1) / numPrepareWorkers;

	std::deque<std::thread> prepareThreads;
	for (int i = 0; i < numPrepareWorkers; ++i) {
		int workerStart = startIndex + workChunks * i;
		int workerEnd = std::min(startIndex + workChunks * (i + 1), endIndex);
		prepareThreads.push_back(std::thread(ipu_prepare, std::cref(swconfig), std::cref(ipuconfig), workerStart, workerEnd, std::cref(A), std::cref(B), std::ref(bufferCmps), std::ref(input_bufs), std::ref(mapping_bufs)));
	}
	
	while (true) {
		int c;
		bufferSelectionMutex.lock();
		for (c = 0; c < bufferCmps.size(); ++c) {
			if (bufferCmps[c] > 0) break;
		}
		bufferSelectionMutex.unlock();
		if (c == bufferCmps.size()) {
			if (startedWorkers == 0 || startedWorkers > stoppedWorkers) {
				continue;
			} else {
				break;
			}
		}
		int numCmps = bufferCmps[c];
		PLOGW << "Getting " << numCmps << " from " << c;
		t.tick();
    auto slot = driver.queue_slot(bufferMaxBucket[c]);
		if (driver.use_remote_buffer) {
			driver.upload(&*input_bufs[c].begin(), &*input_bufs[c].end(), slot);
		}
		driver.prepared_remote_compare(&*input_bufs[c].begin(), &*input_bufs[c].end(), &*results_buf.begin(), &*results_buf.end(), slot);
		t.tock();
		driver.transferResults(
			&*results_buf.begin(), &*results_buf.end(),
			&*mapping_bufs[c].begin(), &*mapping_bufs[c].end(),
			&*(scores.begin() + results_offset), &*(scores.begin() + results_offset + numCmps),
			&*(arange.begin() + results_offset), &*(arange.begin() + results_offset + numCmps),
			&*(brange.begin() + results_offset), &*(brange.begin() + results_offset + numCmps), numCmps);
		results_offset += numCmps;
		bufferSelectionMutex.lock();
		bufferCmps[c] = 0;
		processedBatches += 1;
		bufferSelectionMutex.unlock();
	}

	PLOGW << "Processed " << results_offset << " comparison from " << A.size() << " submitted";

	for (int i = 0; i < prepareThreads.size(); ++i) {
		prepareThreads[i].join();
	}
}

void run_comparison(IpuSwConfig config, std::string referencePath, std::string queryPath) {
	swatlib::TickTock outer; 
	std::vector<swatlib::TickTock> inner(config.numDevices);
	std::vector<ipu::batchaffine::SWAlgorithm> drivers;
	std::deque<std::thread> creatorThreads;
	for (int n = 0; n < config.numDevices; ++n) {
		creatorThreads.push_back(std::thread(ipu_init, std::cref(config), std::ref(drivers)));
	}

	std::vector<std::string> references, queries;
	load_data(referencePath, references);
	load_data(queryPath, queries);

	for (int n = 0; n < config.numDevices; ++n) {
		creatorThreads[n].join();
	}

	PLOGI << "Starting comparisons";

	std::deque<std::thread> ipuThreads;
	outer.tick();
	for (int n = 0; n < config.numDevices; ++n) {
		ipuThreads.push_back(std::thread(ipu_run, std::ref(drivers[n]), std::cref(references), std::cref(queries), n, config.numDevices, std::ref(inner[n])));
	}

	double innerTime = 0;
	for (int n = 0; n < config.numDevices; ++n) {
		ipuThreads[n].join();
		innerTime += static_cast<double>(inner[n].accumulate_microseconds()) / 1e6;
	}
	outer.tock();

	PLOGW << "Submitted batches " << submittedBatches << " processed " << processedBatches;

	uint64_t cellCount = 0;
	for (int i = 0; i < references.size(); ++i) {
		cellCount += references[i].size() * queries[i].size();
	}

	double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
	PLOGI << "outer: " << outerTime << "s inner: " << innerTime << "s";
	double gcupsOuter = static_cast<double>(cellCount) / outerTime / 1e9;
	double gcupsInner = static_cast<double>(cellCount) / innerTime / 1e9;
	PLOGI << "GCUPS outer: " << gcupsOuter << " inner: " << gcupsInner;
}

#endif