#include <iostream>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "cpuswconfig.hpp"

#include "ipuma.h"
#include "cmd_arguments.hpp"
#include "alignment_seqan.hpp"
#include "alignment_genometools.hpp"

using json = nlohmann::json;

template<typename C>
std::vector<int> runAlignment(const ipu::RawSequences& seqs, const ipu::Comparisons& cmps, const ipu::SWConfig& config) {
  swatlib::TickTock t;
  t.tick();
  double cells = 0;

  std::vector<int> scores(cmps.size());

  #pragma omp parallel for
  for (int i = 0; i < cmps.size(); ++i) {
    const auto& cmp = cmps[i];
    cells += (seqs[cmp.indexA].size() * seqs[cmp.indexB].size()) / 1e9;
    // PLOGE << json{
    //   {"i", i},
    //   {"lenH", seqs[cmp.indexA].size()},
    //   {"lenV", seqs[cmp.indexB].size()},
    //   {"seedH", cmp.seedAStartPos},
    //   {"seedV", cmp.seedBStartPos},
    // }.dump();
    int maxScore = 0;
    for (int j = 0; j < NSEEDS; ++j) {
      maxScore = std::max(
				C::align(seqs[cmp.indexA], seqs[cmp.indexB], cmp.seeds[j].seedAStartPos, cmp.seeds[j].seedBStartPos, config.seedLength, config),
				maxScore
			);
    }
    scores[i] = maxScore;
  }
  t.tock();
  double gcups = cells / t.accumulate_microseconds() * 1e6;
  PLOGI << "GCUPS " << gcups;
  return scores;
}

int main(int argc, char** argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);

	cxxopts::Options options("cpusw", "CPU Xdrop implementation");

	options.add_options()
		("c,config", "Configuration file.", cxxopts::value<std::string>())
		("h,help", "Print usage")
		;

	json configJson = CpuSwConfig();
	addArguments(configJson, options, "");

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << "\n";
		exit(0);
	}

	if (result.count("config")) {
		std::string configPath = result["config"].as<std::string>();
		std::ifstream cf(configPath);
		json cj;
		cf >> cj;
		configJson = cj;
	}

	parseArguments(configJson, result);

	CpuSwConfig config = configJson.get<CpuSwConfig>();

	PLOGI << "CPUSWCONFIG" << json{config}.dump();
	auto seqdb = config.getSequences();
	auto [seqs, cmps] = seqdb->get();
	PLOGI << ipu::getDatasetStats(seqs, cmps).dump();

	std::vector<int> scores;
	switch (config.algoconfig.algo) {
	case cpu::Algo::seqan:
		scores = runAlignment<cpu::SeqanAligner>(seqs, cmps, config.swconfig);
		break;
	case cpu::Algo::genometools:
		scores = runAlignment<cpu::GenomeToolsAligner>(seqs, cmps, config.swconfig);
		break;
	}
	if (config.output != "") {
		std::ofstream ofile(config.output);
		ofile << json{scores};
	}

	return 0;

}