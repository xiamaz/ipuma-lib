#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include "ipuswconfig.hpp"
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <omp.h>

#include <plog/Log.h>

#define IPU_JSON_LOG_TAG "IPUSWLOG"

using json = nlohmann::json;

using Fastas = std::vector<swatlib::Fasta>;

std::atomic<int> totalCmpsProcessed = 0;
std::atomic<int> totalBatchesProcessed = 0;

class Barrier {
public:
    explicit Barrier(std::size_t iCount) : 
      mThreshold(iCount), 
      mCount(iCount), 
      mGeneration(0) {
    }

    void wait()
    {
			  std::unique_lock<std::mutex> lLock{mMutex};
        auto lGen = mGeneration;
        if (!--mCount) {
            mGeneration++;
            mCount = mThreshold;
            mCond.notify_all();
        } else {
            mCond.wait(lLock, [this, lGen] { return lGen != mGeneration; });
        }
    }

 private:
     std::mutex mMutex;
    std::condition_variable mCond;
    std::size_t mThreshold;
    std::size_t mCount;
    std::size_t mGeneration;
};

void loadSequences(const std::string& path, std::vector<std::string>& sequences) {
  std::ifstream seqFile(path);
  std::string line;
  while (std::getline(seqFile, line)) {
    sequences.push_back(line);
  }
}

void load_data(const std::string& path, std::vector<std::string>& seqs, int count = 0) {
	if (std::equal(path.end() - 4, path.end(), ".txt")) {
		loadSequences(path, seqs);
	} else {
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
}

void ipu_run(const IpuSwConfig& config, const ipu::RawSequences& A, const ipu::RawSequences& B, const int workerId, const int numWorkers, swatlib::TickTock& t, Barrier& barrier, swatlib::TickTock& outer) {
	const auto& ipuconfig = config.ipuconfig;
	const auto& swconfig = config.swconfig;
	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, workerId);
	PLOGW << "Finished attaching to device";

  std::vector<int> mapping_bufs(ipuconfig.getTotalNumberOfComparisons(), 0);
  std::vector<int32_t> input_bufs(ipuconfig.getInputBufferSize32b());
  size_t results_size = ipuconfig.getTotalNumberOfComparisons() * 3;
  std::vector<int32_t> results_buf(results_size);

	std::vector<int32_t> scores(A.size());
	std::vector<int32_t> arange(A.size());
	std::vector<int32_t> brange(A.size());

	int results_offset = 0;

	int seqsPerWorker = (A.size() + numWorkers - 1) / numWorkers;
	int startIndex = seqsPerWorker * workerId;
	int endIndex = std::min(seqsPerWorker * (workerId + 1), static_cast<int>(A.size()));

	int workRange = endIndex - startIndex;

	PLOGW.printf("%d: Processing %d seqs of total %d seqs", workerId, workRange, A.size());

	int numCmps = 0;
	int totalSize = 0;

	ipu::RawSequences::const_iterator aBegin, bBegin;

	int batchCmpLimit = ipuconfig.getTotalNumberOfComparisons() - ipuconfig.maxBatches;
	int slack = ipuconfig.getTotalBufsize32b() * 4 / 10; // only fill buffer up to 90%, hack to make sure we don't fail
	int batchDataLimit = (ipuconfig.getTotalBufsize32b() * 4) - slack;

	aBegin = A.begin() + startIndex;
	bBegin = B.begin() + startIndex;

	auto submitBatch = [&]() {
		auto aEnd = aBegin + numCmps;
		auto bEnd = bBegin + numCmps;

		std::vector<std::string> aBatch(aBegin, aEnd);
		std::vector<std::string> bBatch(bBegin, bEnd);
		auto maxBucket = ipu::batchaffine::SWAlgorithm::prepare_remote(swconfig, ipuconfig, aBatch, bBatch, &*input_bufs.begin(), &*input_bufs.end(), mapping_bufs.data());

		t.tick();
    auto slot = driver.queue_slot(maxBucket);
		if (driver.algoconfig.useRemoteBuffer) {
			driver.upload(&*input_bufs.begin(), &*input_bufs.end(), slot);
		}
		driver.prepared_remote_compare(&*input_bufs.begin(), &*input_bufs.end(), &*results_buf.begin(), &*results_buf.end(), slot);
		t.tock();
		driver.transferResults(
			&*results_buf.begin(), &*results_buf.end(),
			&*mapping_bufs.begin(), &*mapping_bufs.end(),
			&*(scores.begin() + results_offset), &*(scores.begin() + results_offset + numCmps),
			&*(arange.begin() + results_offset), &*(arange.begin() + results_offset + numCmps),
			&*(brange.begin() + results_offset), &*(brange.begin() + results_offset + numCmps), numCmps);
		results_offset += numCmps;
		totalBatchesProcessed++;
		totalCmpsProcessed += numCmps;
	};

	barrier.wait();
	if (workerId) outer.tick();
	for (int i = startIndex; i < endIndex; ++i) {
		const auto alen = A[i].size();
		const auto blen = B[i].size();
		if ((numCmps + 1 < batchCmpLimit) && (totalSize + alen + blen < batchDataLimit)) {
			numCmps++;
			totalSize += alen + blen;
		} else {
			submitBatch();

			aBegin = aBegin + numCmps + 1;
			bBegin = bBegin + numCmps + 1;
			numCmps = 1;
			totalSize = alen + blen;
		}
	}
	if (numCmps > 0) {
		submitBatch();
	}
	barrier.wait();
	if (workerId) outer.tock();
}

void run_comparison(IpuSwConfig config, std::string referencePath, std::string queryPath) {
	std::vector<std::string> references, queries;
	swatlib::TickTock outer; 
	std::vector<swatlib::TickTock> inner(config.numDevices);
	// std::vector<ipu::batchaffine::SWAlgorithm> drivers;
	std::deque<std::thread> creatorThreads;

	json configLog = {
		{"tag", "run_comparison_setup"},
		{"config", config},
		{"ref_path", referencePath},
		{"query_path", queryPath},
	};
	PLOGW << IPU_JSON_LOG_TAG << configLog.dump();

	load_data(referencePath, references);
	load_data(queryPath, queries);

	if (references.size() != queries.size()) {
		PLOGW << "refsize " << references.size();
		PLOGW << "quersize " << queries.size();
		exit(1);
	}

	PLOGI << "Starting comparisons";

	Barrier barrier(config.numDevices);

	std::deque<std::thread> ipuThreads;
	for (int n = 0; n < config.numDevices; ++n) {
		ipuThreads.push_back(std::thread(ipu_run, std::cref(config), std::cref(references), std::cref(queries), n, config.numDevices, std::ref(inner[n]), std::ref(barrier), std::ref(outer)));
	}


	double innerTime = 0;
	double maxInnerTime = 0;
	for (int n = 0; n < config.numDevices; ++n) {
		ipuThreads[n].join();
		auto threadSeconds = static_cast<double>(inner[n].accumulate_microseconds()) / 1e6;
		innerTime += threadSeconds;
		maxInnerTime = std::max(maxInnerTime, threadSeconds);
	}

	PLOGW << "Total cmps processed/submitted " << totalCmpsProcessed << "/" << references.size();
	PLOGW << "Total Batches Processed " << totalBatchesProcessed;

	uint64_t cellCount = 0;
	for (int i = 0; i < references.size(); ++i) {
		cellCount += references[i].size() * queries[i].size();
	}

	double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
	double gcupsOuter = static_cast<double>(cellCount) / outerTime / 1e9;
	double gcupsInner = static_cast<double>(cellCount) / innerTime / 1e9;
	double gcupsInnerMaxTime = static_cast<double>(cellCount) / maxInnerTime / 1e9 / config.numDevices;
	json logdata = {
		{"tag", "final_log"},
		{"cells", cellCount},
		{"outer_time_s", outerTime},
		{"inner_time_acc_s", innerTime},
		{"inner_time_max_s", maxInnerTime},
		{"outer_gcups", gcupsOuter},
		{"inner_gcups_acc", gcupsInner},
		{"inner_gcups_max_time", gcupsInnerMaxTime},
	};
	PLOGW << IPU_JSON_LOG_TAG << logdata.dump();
}

#endif