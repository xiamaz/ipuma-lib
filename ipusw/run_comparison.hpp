#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include "ipuswconfig.hpp"
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>

using json = nlohmann::json;

using Fastas = std::vector<swatlib::Fasta>;

void load_data(const std::string& path, Fastas& seqs, int count = 0) {
  std::ifstream is;
  is.open(path);
  if (is.fail()) {
      throw std::runtime_error("Opening file at " + path + " failed.");
  }
  swatlib::Fasta f;
	int i = 0;
	while (!(count) || i < count) {
      if (is >> f) {
          seqs.push_back(f);
      } else {
          seqs.push_back(f);
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

	Fastas referenceFasta, queryFasta;
	PLOGI << "Loading data from " << referencePath;
	load_data(referencePath, referenceFasta);
	PLOGI << "Loading data from " << queryPath;
	load_data(queryPath, queryFasta);

	std::vector<std::string> references, queries;
	ipu::partition::BucketMap map(config.ipuconfig.tilesUsed, config.ipuconfig.maxBatches, config.ipuconfig.bufsize);
	int batchCmps = 0, curBucket = 0;
  ipu::partition::BucketHeap q;
	if (config.ipuconfig.fillAlgo == ipu::Algorithm::greedy) {
  	for (auto& b : map.buckets) {
  	  q.push(std::ref(b));
  	}
	}

	PLOGI << "Starting comparisons";

	std::vector<int32_t> scores(referenceFasta.size());
	std::vector<int32_t> arange(referenceFasta.size());
	std::vector<int32_t> brange(referenceFasta.size());

	outer.tick();
	int results_offset = 0;
	for (int i = 0; i < referenceFasta.size(); ++i) {
		const auto& a = referenceFasta[i].sequence;
		const auto& b = queryFasta[i].sequence;
		if (ipu::partition::fillBuckets(config.ipuconfig.fillAlgo, map, {a}, {b}, batchCmps, curBucket, q)) {
			batchCmps++;
			references.push_back(a);
			queries.push_back(b);
		} else {
			driver.prepare_remote(config.swconfig, config.ipuconfig, references, queries, &*inputs_buf.begin(), &*inputs_buf.end(), mapping.data());

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

			results_offset += batchCmps;
			references = {a};
			queries = {b};
			batchCmps = 1;
			curBucket = 0;
			map = ipu::partition::BucketMap(config.ipuconfig.tilesUsed, config.ipuconfig.maxBatches, config.ipuconfig.bufsize);
			q = {};
			if (config.ipuconfig.fillAlgo == ipu::Algorithm::greedy) {
  			for (auto& b : map.buckets) {
  			  q.push(std::ref(b));
  			}
			}
			if (!ipu::partition::fillBuckets(config.ipuconfig.fillAlgo, map, {a}, {b}, batchCmps, curBucket, q)) {
				PLOG_FATAL << "Really out of buckets.\n";
			}
		}
	}
	if (batchCmps > 0) {
		driver.prepare_remote(config.swconfig, config.ipuconfig, references, queries, &*inputs_buf.begin(), &*inputs_buf.end(), mapping.data());

		inner.tick();
		if (driver.use_remote_buffer) {
			ipu::batchaffine::slotToken slot = driver.upload(&*inputs_buf.begin(), &*inputs_buf.end());
		}
		driver.prepared_remote_compare(&*inputs_buf.begin(), &*inputs_buf.end(), &*results_buf.begin(), &*results_buf.end());
		inner.tock();
	}
	outer.tock();

	uint64_t cellCount = 0;
	for (int i = 0; i < referenceFasta.size(); ++i) {
		cellCount += referenceFasta[i].sequence.size() * queryFasta[i].sequence.size();
	}

	double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
	double innerTime = static_cast<double>(inner.accumulate_microseconds()) / 1e6;
	PLOGI << "outer: " << outerTime << "s inner: " << innerTime << "s";
	double gcupsOuter = static_cast<double>(cellCount) / outerTime / 1e9;
	double gcupsInner = static_cast<double>(cellCount) / innerTime / 1e9;
	PLOGI << "GCUPS outer: " << gcupsOuter << " inner: " << gcupsInner;
}

#endif