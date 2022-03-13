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

void loadSequences(const std::string &path, std::vector<std::string> &sequences)
{
  std::ifstream seqFile(path);
  std::string line;
  while (std::getline(seqFile, line))
  {
    sequences.push_back(line);
  }
}

void load_data(const std::string &path, std::vector<std::string> &seqs, int count = 0)
{
  if (std::equal(path.end() - 4, path.end(), ".txt"))
  {
    PLOGI << "Loading all entries from " << path;
    loadSequences(path, seqs);
  }
  else
  {
    if (count > 0)
    {
      PLOGI << "Loading " << count << " entries from " << path;
    }
    else
    {
      PLOGI << "Loading all entries from " << path;
    }
    std::ifstream is;
    is.open(path);
    if (is.fail())
    {
      throw std::runtime_error("Opening file at " + path + " failed.");
    }
    swatlib::Fasta f;
    int i = 0;
    while (!(count) || i < count)
    {
      if (is >> f)
      {
        seqs.push_back(f.sequence);
      }
      else
      {
        seqs.push_back(f.sequence);
        break;
      }
      ++i;
    }
  }
}

void run_comparison(IpuSwConfig config, std::string referencePath, std::string queryPath) {
  std::vector<std::string> references, queries;
  swatlib::TickTock outer;
  std::vector<swatlib::TickTock> inner(config.numThreads);

  int duplicationFactor = config.duplicateDatasets ? std::max(config.numDevices / 2, 1) : 1;

  json configLog = {
      {"tag", "run_comparison_setup"},
      {"config", config},
      {"ref_path", referencePath},
      {"query_path", queryPath},
      {"duplicationFactor", duplicationFactor}
  };
  PLOGW << IPU_JSON_LOG_TAG << configLog.dump();

  load_data(referencePath, references);
  load_data(queryPath, queries);

  const int batchCmpLimit = config.ipuconfig.getTotalNumberOfComparisons() - config.ipuconfig.maxBatches;
  const int batchDataLimit = config.ipuconfig.getTotalBufsize32b() * 4 - config.ipuconfig.bufsize * 100;
  auto batches = createBatches(references, queries, batchCmpLimit, batchDataLimit);

  if (config.duplicateDatasets) {
    int originalSize = batches.size();
    int duplicatedSize = duplicationFactor * originalSize;
    batches.resize(duplicatedSize);
    for (int i = originalSize; i < duplicatedSize; ++i) {
      // copy older batch data into new data
      batches[i] = batches[i % originalSize];
    }
    PLOGW.printf("Duplicating dataset %dx from %d to %d", duplicationFactor, originalSize, duplicatedSize);
  }

  auto driver = IPUMultiDriver(config);
  std::vector<std::thread> workerThreads;
  PLOGI << "Starting comparisons";
  outer.tick();
  for (int n = 0; n < config.numThreads; ++n) {
    workerThreads.push_back(std::thread(runIpuWorker, n, config.numThreads, std::ref(driver), std::cref(references), std::cref(queries), std::cref(batches), std::ref(inner[n])));
  }

  double innerTime = 0;
  double maxInnerTime = 0;
  for (int n = 0; n < config.numThreads; ++n) {
    workerThreads[n].join();
    auto threadSeconds = static_cast<double>(inner[n].accumulate_microseconds()) / 1e6;
    innerTime += threadSeconds;
  }
  outer.tock();

  PLOGW << "Total cmps processed/submitted " << driver.totalCmpsProcessed << "/" << references.size();
  PLOGW << "Total Batches Processed " << driver.totalBatchesProcessed;

  uint64_t cellCount = 0;
  for (int i = 0; i < references.size(); ++i)
  {
    cellCount += references[i].size() * queries[i].size();
  }
  
  if (config.duplicateDatasets) {
    // cave: we are only duplicating the batches currently. This needs to be removed if the input data itself is duplicated
    cellCount *= duplicationFactor;
  }

  double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
  double gcupsOuter = static_cast<double>(cellCount) / outerTime / 1e9;
  double gcupsInner = static_cast<double>(cellCount) / innerTime / 1e9;
  json logdata = {
      {"tag", "final_log"},
      {"cells", cellCount},
      {"outer_time_s", outerTime},
      {"inner_time_acc_s", innerTime},
      {"inner_time_max_s", maxInnerTime},
      {"outer_gcups", gcupsOuter},
      {"inner_gcups_acc", gcupsInner},
  };
  PLOGW << IPU_JSON_LOG_TAG << logdata.dump();
}

#endif