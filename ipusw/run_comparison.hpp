#ifndef RUN_COMPARISON_HPP
#define RUN_COMPARISON_HPP
#include "ipuswconfig.hpp"
#include <string>

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
	Fastas referenceFasta, queryFasta;
	load_data(referencePath, referenceFasta);
	load_data(queryPath, queryFasta);

	std::vector<std::string> references, queries;
	for (int i = 0; i < referenceFasta.size(); ++i) {
		references.push_back(referenceFasta[i].sequence);
	}
	for (int i = 0; i < queryFasta.size(); ++i) {
		queries.push_back(queryFasta[i].sequence);
	}

	auto driver = ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig);
	driver.compare_local(references, queries);
}

#endif