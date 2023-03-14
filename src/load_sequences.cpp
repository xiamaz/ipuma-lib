#include<tuple>
#include<istream>
#include<fstream>
#include <nlohmann/json.hpp>

#include <plog/Log.h>

#include "types.h"
#include "load_sequences.h"

using json = nlohmann::json;
namespace ipu {
json getDatasetStats(const ipu::RawSequences& seqs, const ipu::Comparisons& cmps) {
  auto numComparisons = cmps.size();
  size_t maxSequenceLength = 0;
  for (const auto& s : seqs) {
    maxSequenceLength = std::max(maxSequenceLength, s.size());
  }

  json stats = {
    {"numComparisons", numComparisons},
    {"maxSequenceLength", maxSequenceLength},
  };
  return stats;

  PLOGI << "COMPARISON STATS: " << stats.dump();
}

std::unique_ptr<std::ifstream> openFile(std::string path) {
  std::unique_ptr<std::ifstream> file(new std::ifstream(path));
  if (!file->good()) {
    PLOGE << "Failed to open " << path;
    file->close();
    return {nullptr};
  }
  PLOGD << "Loading from " << path;
  return file;
}

typedef bool FailFlag;
typedef bool EOFFlag;

std::tuple<std::string, int, int, FailFlag, EOFFlag> readFile(std::unique_ptr<std::ifstream>& seqF, std::unique_ptr<std::ifstream>& seed1F, std::unique_ptr<std::ifstream>& seed2F) {
  std::string s = "";
  int seed1 = -1, seed2 = -1;
  std::getline(*seqF, s);
  bool allEof = seqF->eof();
  bool anyFail = seqF->fail();
  if (seed1F) {
    *seed1F >> seed1;
    allEof = allEof && seed1F->eof();
    anyFail = anyFail || seed1F->fail();
  }
  if (seed2F) {
    *seed2F >> seed2;
    allEof = allEof && seed2F->eof();
    anyFail = anyFail || seed2F->fail();
  }

  return {std::move(s), seed1, seed2, anyFail, allEof};
}

SequenceData::SequenceData(const SequenceConfig& config) {
  auto fileH = openFile(config.seqsH);
  auto fileV = openFile(config.seqsV);
  auto fileHs1 = openFile(config.seedsH1);
  auto fileVs1 = openFile(config.seedsV1);
  auto fileHs2 = openFile(config.seedsH2);
  auto fileVs2 = openFile(config.seedsV2);
  int sH, sV;
  std::string lH, lV;
  int i = 0;
  while (true) {
    auto [lH, sH1, sH2, anyFailH, eofH] = readFile(fileH, fileHs1, fileHs2);
    auto [lV, sV1, sV2, anyFailV, eofV] = readFile(fileV, fileVs1, fileVs2);
    if (eofH && eofV) {
      break;
    }
    if (!anyFailH && !anyFailV) {
      sequences.push_back(lH);
      seqs.push_back(sequences.back());
      sequences.push_back(lV);
      seqs.push_back(sequences.back());
      const int sizeA = (int) seqs[(int) (2 * i)].size();
      const int sizeB = (int) seqs[(int) (2 * i + 1)].size();
      cmps.push_back({
        i,
        (int) (2 * i),
        sizeA,
        (int) (2 * i + 1),
        sizeB,
        {{{sH1, sH2}, {sV1, sV2}}},
      });
      i++;
    } else {
      if (anyFailH) {
        PLOGE << "Malformed sequence input in H sequence string or seeds";
      }
      if (anyFailV) {
        PLOGE << "Malformed sequence input in V sequence string or seeds";
      }
      PLOGF << "Failed to read input files";
      exit(1);
    }
  }
}

void from_json(const json& j, SequenceConfig& c) {
	j.at("seqsH").get_to(c.seqsH);
	j.at("seqsV").get_to(c.seqsV);
	j.at("seedsH1").get_to(c.seedsH1);
	j.at("seedsH2").get_to(c.seedsH2);
	j.at("seedsV1").get_to(c.seedsV1);
	j.at("seedsV2").get_to(c.seedsV2);
	j.at("complexity1").get_to(c.complexity1);
	j.at("complexity2").get_to(c.complexity2);
}

void to_json(json& j, const SequenceConfig& c) {
	j = json{
		{"seqsH", c.seqsH},
		{"seqsV", c.seqsV},
		{"seedsH1", c.seedsH1},
		{"seedsH2", c.seedsH2},
		{"seedsV1", c.seedsV1},
		{"seedsV2", c.seedsV2},
		{"complexity1", c.complexity1},
		{"complexity2", c.complexity2},
	};
}


std::tuple<ipu::RawSequences, ipu::Comparisons> SequenceData::get() {
  return {std::move(seqs), std::move(cmps)};
}
}