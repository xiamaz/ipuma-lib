#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include "ipuswconfig.hpp"
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <plog/Log.h>

using json = nlohmann::json;

using Fastas = std::vector<swatlib::Fasta>;

void load_data(const std::string& path, std::vector<std::string>& seqs, int count = 0) {
	PLOGI << "Loading " << count << " entries from " << path;
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

void run_comparison(IpuSwConfig config, std::string referencePath, std::string queryPath) {
  std::vector<int> mapping(config.ipuconfig.getTotalNumberOfComparisons(), 0);
  size_t inputs_size = config.ipuconfig.getInputBufferSize32b();
  std::vector<int32_t> inputs_buf(inputs_size);
  size_t results_size = config.ipuconfig.getTotalNumberOfComparisons() * 3;
  std::vector<int32_t> results_buf(results_size);

	swatlib::TickTock outer, inner;
	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, false);

	std::vector<std::string> references, queries;
	load_data(referencePath, references);
	load_data(queryPath, queries);

	ipu::partition::BucketMap map(config.ipuconfig.tilesUsed, config.ipuconfig.maxBatches, config.ipuconfig.bufsize);
	int batchCmps = 0, curBucket = 0;
  ipu::partition::BucketHeap q;
	if (config.ipuconfig.fillAlgo == ipu::Algorithm::greedy) {
  	for (auto& b : map.buckets) {
  	  q.push(std::ref(b));
  	}
	}

	PLOGI << "Starting comparisons";

	std::vector<int32_t> scores(references.size());
	std::vector<int32_t> arange(references.size());
	std::vector<int32_t> brange(references.size());

	outer.tick();
	int results_offset = 0;

	// heuritic slicing for greater speed
	int numCmps = 0;
	int totalSize = 0;

	int batchCmpLimit = config.ipuconfig.getTotalNumberOfComparisons() - config.ipuconfig.maxBatches;
	int batchDataLimit = (config.ipuconfig.getTotalBufsize32b() * 4) - config.ipuconfig.bufsize;

	ipu::RawSequences::iterator aBegin, aEnd, bBegin, bEnd;
	aBegin = references.begin();
	bBegin = queries.begin();

	for (int i = 0; i < references.size(); ++i) {
		const auto alen = references[i].size();
		const auto blen = queries[i].size();
		if (numCmps + 1 < batchCmpLimit && totalSize + alen + blen < batchDataLimit) {
			numCmps++;
			totalSize += alen + blen;
		} else {
			aEnd = aBegin + numCmps;
			bEnd = bBegin + numCmps;

			std::vector<std::string> aBatch(aBegin, aEnd);
			std::vector<std::string> bBatch(bBegin, bEnd);
			driver.prepare_remote(config.swconfig, config.ipuconfig, aBatch, bBatch, &*inputs_buf.begin(), &*inputs_buf.end(), mapping.data());

			inner.tick();
			if (driver.use_remote_buffer) {
				ipu::batchaffine::slotToken slot = driver.upload(&*inputs_buf.begin(), &*inputs_buf.end());
			}
			driver.prepared_remote_compare(&*inputs_buf.begin(), &*inputs_buf.end(), &*results_buf.begin(), &*results_buf.end());
			inner.tock();

			driver.transferResults(
				&*results_buf.begin(), &*results_buf.end(),
				&*mapping.begin(), &*mapping.end(),
				&*(scores.begin() + results_offset), &*(scores.begin() + results_offset + batchCmps),
				&*(arange.begin() + results_offset), &*(arange.begin() + results_offset + batchCmps),
				&*(brange.begin() + results_offset), &*(brange.begin() + results_offset + batchCmps), batchCmps);

			results_offset += numCmps;
			numCmps = 1;
			totalSize = alen + blen;
			aBegin = aEnd + 1;
			bBegin = bEnd + 1;
		}
	}
	if (batchCmps > 0) {
		std::vector<std::string> aBatch(aBegin, aEnd);
		std::vector<std::string> bBatch(bBegin, bEnd);
		driver.prepare_remote(config.swconfig, config.ipuconfig, aBatch, bBatch, &*inputs_buf.begin(), &*inputs_buf.end(), mapping.data());

		inner.tick();
		if (driver.use_remote_buffer) {
			ipu::batchaffine::slotToken slot = driver.upload(&*inputs_buf.begin(), &*inputs_buf.end());
		}
		driver.prepared_remote_compare(&*inputs_buf.begin(), &*inputs_buf.end(), &*results_buf.begin(), &*results_buf.end());
		inner.tock();
		driver.transferResults(
			&*results_buf.begin(), &*results_buf.end(),
			&*mapping.begin(), &*mapping.end(),
			&*(scores.begin() + results_offset), &*(scores.begin() + results_offset + batchCmps),
			&*(arange.begin() + results_offset), &*(arange.begin() + results_offset + batchCmps),
			&*(brange.begin() + results_offset), &*(brange.begin() + results_offset + batchCmps), batchCmps);
	}
	outer.tock();

	uint64_t cellCount = 0;
	for (int i = 0; i < references.size(); ++i) {
		cellCount += references[i].size() * queries[i].size();
	}

	double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
	double innerTime = static_cast<double>(inner.accumulate_microseconds()) / 1e6;
	PLOGI << "outer: " << outerTime << "s inner: " << innerTime << "s";
	double gcupsOuter = static_cast<double>(cellCount) / outerTime / 1e9;
	double gcupsInner = static_cast<double>(cellCount) / innerTime / 1e9;
	PLOGI << "GCUPS outer: " << gcupsOuter << " inner: " << gcupsInner;
}

#endif