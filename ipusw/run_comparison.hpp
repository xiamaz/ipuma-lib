#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include "ipuswconfig.hpp"
#include "ipumultidriver.hpp"
#include "ipuworker.hpp"
#include "batch.hpp"
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

void loadSequences(const std::string &path, std::vector<std::string> &sequences) {
  std::ifstream seqFile(path);
  if (!seqFile.good()) {
    PLOGF << "File not found: " << path;
    exit(1);
  }
  std::string line;
  while (std::getline(seqFile, line)) {
    sequences.push_back(line);
  }
}

void load_data(const std::string &path, std::vector<std::string> &seqs, int count = 0) {
  if (std::equal(path.end() - 4, path.end(), ".txt")) {
    PLOGI << "Loading all entries from " << path;
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

void init_driver(IPUMultiDriver& driver) {
  driver.init();
}

void run_comparison(IpuSwConfig config, std::string referencePath, std::string queryPath) {
  auto driver = IPUMultiDriver(config);
  std::vector<std::string> references, queries;
  swatlib::TickTock outer;

  int duplicationFactor = config.duplicateDatasets ? std::max(config.numDevices, 1) : 1;
  if (config.duplicationFactor > 0) {
    duplicationFactor = config.duplicationFactor;
  }

  json configLog = {
      {"tag", "run_comparison_setup"},
      {"config", config},
      {"ref_path", referencePath},
      {"query_path", queryPath},
      {"duplicationFactor", duplicationFactor}
  };
  PLOGW << IPU_JSON_LOG_TAG << configLog.dump();

  std::thread refLoader(load_data, std::cref(referencePath), std::ref(references), 0);
  std::thread queryLoader(load_data, std::cref(queryPath), std::ref(queries), 0);
  std::thread driverInit(init_driver, std::ref(driver));

  refLoader.join();
  queryLoader.join();

  const int batchCmpLimit = config.ipuconfig.getTotalNumberOfComparisons() - config.ipuconfig.maxBatches;
  const int batchDataLimit = config.ipuconfig.getTotalBufsize32b() * 4 - config.ipuconfig.bufsize * 100;
  PLOGF << references.size();
  auto batches = createBatches(references, queries, batchCmpLimit, batchDataLimit);

  if (config.duplicateDatasets) {
    int originalSize = batches.size();
    int duplicatedSize = duplicationFactor * originalSize;
    PLOGW.printf("Duplicating dataset %dx from %d to %d", duplicationFactor, originalSize, duplicatedSize);
  }
  if (batches.size() > WORK_QUEUE_SIZE) {
    PLOGF.printf("Batches are larger than work queue size: %d > %d", batches.size(), WORK_QUEUE_SIZE);
  }

  std::vector<std::thread> workerThreads;
  for (int n = 0; n < config.numThreads; ++n) {
    workerThreads.push_back(std::thread(runIpuWorker, n, config.numThreads, std::ref(driver), std::cref(references), std::cref(queries), std::ref(batches)));
  }
  for (int n = 0; n < config.numThreads; ++n) {
    workerThreads[n].join();
  }

  std::vector<BatchResult> results(batches.size() * duplicationFactor);
	const auto resultBufferSize = driver.config.ipuconfig.getTotalNumberOfComparisons() * 3;
  for (auto& r : results) {
    r.resultBuffer.resize(resultBufferSize);
  }

  driverInit.join();
  driver.run();

  PLOGI << "Starting comparisons";
  std::vector<std::thread> receiverThreads;
  const int receiverThreadNum = config.numThreads / 2;
  // std::vector<ipu::batchaffine::Job*> jobs(batches.size() * duplicationFactor, nullptr);
  for (int i = 0; i < receiverThreadNum; ++i) {
    receiverThreads.push_back(std::thread(workerResult, i, receiverThreadNum, std::ref(driver), std::ref(batches), std::ref(results)));
  }
  outer.tick();
  for (int i = 0; i < batches.size() * duplicationFactor; ++i) {
    auto& batch = batches[i % batches.size()];
    auto& result = results[i];
    // batch.job = driver.submit(batch.inputBuffer, batch.resultBuffer);
    auto* job = driver.submit(batch.inputBuffer, result.resultBuffer);
    results[i].job = job;
    // jobs[i] = job;
    driver.totalCmpsSubmitted += batch.numCmps;
    driver.totalBatchesSubmitted += 1;
  }
  for (int i = 0; i < receiverThreadNum; ++i) {
    receiverThreads[i].join();
  }
  outer.tock();

  uint64_t cellCount = 0;
  for (int i = 0; i < references.size(); ++i) {
    cellCount += references[i].size() * queries[i].size();
  }
  
  if (config.duplicateDatasets) {
    // cave: we are only duplicating the batches currently. This needs to be removed if the input data itself is duplicated
    cellCount *= duplicationFactor;
  }

  double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
  double gcupsOuter = static_cast<double>(cellCount) / outerTime / 1e9;
  json logdata = {
      {"tag", "final_log"},
      {"cells", cellCount},
      {"outer_time_s", outerTime},
      {"outer_gcups", gcupsOuter},
      {"cmps_processed", driver.totalCmpsProcessed.load()},
      {"cmps_submitted", driver.totalCmpsSubmitted.load()},
      {"cmps_orig", references.size()},
      {"batches_processed", driver.totalBatchesProcessed.load()},
      {"batches_submitted", driver.totalBatchesSubmitted.load()},
      {"batches_orig", batches.size()},
  };
  PLOGW << IPU_JSON_LOG_TAG << logdata.dump();
}

#endif