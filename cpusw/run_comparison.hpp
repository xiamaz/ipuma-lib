#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include <string>

#include "cpuswconfig.hpp"
#include "aligners.hpp"

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

void run_comparison(CpuSwConfig config, std::string referencePath, std::string queryPath) {
	Fastas referenceFasta, queryFasta;
	PLOGI << "Loading data from " << referencePath;
	load_data(referencePath, referenceFasta);
	PLOGI << "Loading data from " << queryPath;
	load_data(queryPath, queryFasta);

	PLOGI << "Casting into strings.";
	std::vector<std::string> references, queries;
	for (int i = 0; i < referenceFasta.size(); ++i) {
		references.push_back(referenceFasta[i].sequence);
	}
	for (int i = 0; i < queryFasta.size(); ++i) {
		queries.push_back(queryFasta[i].sequence);
	}

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