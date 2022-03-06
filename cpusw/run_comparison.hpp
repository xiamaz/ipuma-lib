#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include <string>

#include "cpuswconfig.hpp"
#include "aligners.hpp"

using json = nlohmann::json;

using Fastas = std::vector<swatlib::Fasta>;

void loadSequences(const std::string& path, std::vector<std::string>& sequences) {
  std::ifstream seqFile(path);
  std::string line;
  while (std::getline(seqFile, line)) {
    sequences.push_back(line);
  }
}

void load_data(const std::string& path, std::vector<std::string>& seqs, int count = 0) {
	if (std::equal(path.end() - 4, path.end(), ".txt")) {
		PLOGI << "Loading all entries from " << path;
		loadSequences(path, seqs);
	} else {
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

void run_comparison(CpuSwConfig config, std::string referencePath, std::string queryPath) {
	std::vector<std::string> references, queries;
	PLOGI << "Loading data from " << referencePath;
	load_data(referencePath, references);
	PLOGI << "Loading data from " << queryPath;
	load_data(queryPath, queries);

	PLOGI << "Compare using " << static_cast<json>(config.algoconfig.algo).dump();
	std::unique_ptr<CPUAligner> al(nullptr);
	switch (config.algoconfig.algo) {
	case (cpu::Algo::stripedsw):
		al.reset(new SSWAligner(config));
		break;
	}
	al->compare(references, queries);

	auto results = al->get_results();
	PLOGI << "Num cmps: " << results.scores.size();
	PLOGI << "Max score: " << *std::max_element(results.scores.begin(), results.scores.end());
}

#endif