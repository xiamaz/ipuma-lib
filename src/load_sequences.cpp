#include<tuple>
#include<istream>
#include<fstream>
#include <nlohmann/json.hpp>

#include <plog/Log.h>

#include "types.h"

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

std::tuple<ipu::RawSequences, ipu::Comparisons> prepareComparisons(std::string seqH, std::string seqV, std::string seedH, std::string seedV) {
  ipu::RawSequences seqs;
  ipu::Comparisons cmps;

  std::ifstream fileH(seqH);
  std::ifstream fileV(seqV);
  std::ifstream fileHs(seedH);
  std::ifstream fileVs(seedV);
  int sH, sV;
  std::string lH, lV;
  int i = 0;
  PLOGD << "Loading from " << seqH << " " << seqV << " " << seedH << " " << seedV;
  if (!fileH.good()) {
    PLOGE << "Failed to open " << seqH;
  }
  if (!fileV.good()) {
    PLOGE << "Failed to open " << seqV;
  }
  if (!fileHs.good()) {
    PLOGE << "Failed to open " << seedH;
  }
  if (!fileVs.good()) {
    PLOGE << "Failed to open " << seedV;
  }
  while (true) {
    std::getline(fileH, lH);
    std::getline(fileV, lV);
    fileHs >> sH;
    fileVs >> sV;
    bool fHValid = !fileH.fail();
    bool fVValid = !fileV.fail();
    bool fHsValid = !fileHs.fail();
    bool fVsValid = !fileVs.fail();
    auto allFilesEOF = fileH.eof() & fileH.eof() & fileHs.eof() & fileVs.eof();
    if (allFilesEOF) {
      break;
    }
    if (fHValid & fVValid & fHsValid & fVsValid) {
      seqs.push_back(lH);
      seqs.push_back(lV);
      cmps.push_back({
        .indexA = (int) (2 * i),
        .indexB = (int) (2 * i + 1),
        .seedAStartPos = (int) sH,
        .seedBStartPos = (int) sV,
      });
      i++;
    } else {
      if (!fHValid) {
        PLOGE << "Malformed sequence input in " << seqH;
      }
      if (!fVValid) {
        PLOGE << "Malformed sequence input in " << seqV;
      }
      if (!fHsValid) {
        PLOGE << "Malformed seed input in " << seedH;
      }
      if (!fVsValid) {
        PLOGE << "Malformed seed input in " << seedV;
      }
      PLOGF << "Failed to read input files";
      exit(1);
    }
  }

  return {std::move(seqs), std::move(cmps)};
}
}