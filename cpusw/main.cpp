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
#include "alignment_libgaba.hpp"
#include "alignment_ksw2.hpp"
#include "alignment_ipuma_cpu.hpp"

using json = nlohmann::json;

template<typename C>
std::vector<int> runAlignment(const ipu::RawSequences& seqs, const ipu::Comparisons& cmps, const ipu::SWConfig& config, int threads) {
  swatlib::TickTock t;
  t.tick();
  double gcells = 0;

  std::vector<int> scores(cmps.size());

	if (threads > 0) {
		omp_set_num_threads(threads);
	}

	C comparator(config);

  #pragma omp parallel for
  for (int i = 0; i < cmps.size(); ++i) {
    const auto& cmp = cmps[i];
    // PLOGE << json{
    //   {"i", i},
    //   {"lenH", seqs[cmp.indexA].size()},
    //   {"lenV", seqs[cmp.indexB].size()},
    //   {"seedH", cmp.seedAStartPos},
    //   {"seedV", cmp.seedBStartPos},
    // }.dump();
    int maxScore = -std::numeric_limits<int>::infinity();
    for (int j = 0; j < NSEEDS; ++j) {
			if ((cmp.seeds[j].seedAStartPos < 0) || (cmp.seeds[j].seedBStartPos < 0)) continue;
    	gcells += (seqs[cmp.indexA].size() * seqs[cmp.indexB].size()) / 1e9;
      maxScore = std::max(
				comparator.align(seqs[cmp.indexA], seqs[cmp.indexB], cmp.seeds[j].seedAStartPos, cmp.seeds[j].seedBStartPos, config.seedLength),
				maxScore
			);
    }
    scores[i] = maxScore;
  }
  t.tock();
	auto time_us = t.accumulate_microseconds();
  double gcups = gcells / time_us * 1e6;

	auto jsonlog = json{
		{"time_ms", time_us / 1e3},
		{"gcups", gcups},
		{"gigacells", gcells},
		{"comparisons", cmps.size()},
	};

	PLOGI << jsonlog.dump();

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
	auto seqdb = config.loaderconfig.getMultiSequences(config.swconfig);
	auto [seqs, mcmps] = seqdb.get();
	PLOGI << ipu::getDatasetStats(seqs, mcmps).dump();

	ipu::Comparisons cmps = convertToComparisons(mcmps);

	std::vector<int> scores;
	switch (config.algoconfig.algo) {
	case cpu::Algo::seqan:
		scores = runAlignment<cpu::SeqanAligner>(seqs, cmps, config.swconfig, config.algoconfig.threads);
		break;
	case cpu::Algo::genometools:
		scores = runAlignment<cpu::GenomeToolsAligner>(seqs, cmps, config.swconfig, config.algoconfig.threads);
		break;
	case cpu::Algo::libgaba:
		scores = runAlignment<cpu::GabaAligner>(seqs, cmps, config.swconfig, config.algoconfig.threads);
		break;
	case cpu::Algo::ksw2:
		scores = runAlignment<cpu::Ksw2Aligner>(seqs, cmps, config.swconfig, config.algoconfig.threads);
		break;
	case cpu::Algo::ipumacpu:
		scores = runAlignment<cpu::IpumaCpuAligner>(seqs, cmps, config.swconfig, config.algoconfig.threads);
		break;
	}
	if (config.output != "") {
		std::ofstream ofile(config.output);
		ofile << json{scores};
	}

	return 0;

}