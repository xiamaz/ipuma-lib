#include<tuple>
#include<istream>
#include<fstream>
#include <nlohmann/json.hpp>

#include <plog/Log.h>

#include "types.h"
#include "sequences_loader.h"

using json = nlohmann::json;
namespace ipu {
template<typename C>
size_t getNumComparisons(const std::vector<C>& cmps);

template<>
size_t getNumComparisons(const std::vector<Comparison>& cmps) {
  return cmps.size();
}

template<>
size_t getNumComparisons(const std::vector<MultiComparison>& cmps) {
  size_t cmpCount = 0;
  for (const auto& c : cmps) {
    cmpCount += c.comparisonCount;
  }
  return cmpCount;
}

template<typename C>
double getGcells(const C& comparison, int seedLen, bool useSeed);

template<>
double getGcells(const Comparison& comparison, int seedLen, bool useSeed) {
  double gcells = 0;
  for (const auto& seed : comparison.seeds) {
    if (seed.seedAStartPos >= 0 && seed.seedBStartPos >= 0) {
      if (useSeed) {
        gcells += (seed.seedAStartPos * seed.seedBStartPos + (
          std::max(comparison.sizeA - (seed.seedAStartPos + seedLen), 0) * std::max(comparison.sizeB - (seed.seedBStartPos + seedLen), 0)
        )) / 1e9;

      } else {
        gcells += comparison.sizeA * comparison.sizeB / 1e9;
      }
    }
  }
  return gcells;
}

template<>
double getGcells(const MultiComparison& comparison, int seedLen, bool useSeed) {
  double gcells = 0;
  for (const auto& c : comparison.comparisons) {
    gcells += getGcells<Comparison>(c, seedLen, useSeed);
  }
  return gcells;
}

template<typename C>
json getDatasetStats(const ipu::RawSequences& seqs, const std::vector<C>& cmps, int seedLen) {
  auto numComparisons = getNumComparisons<C>(cmps);
  size_t maxSequenceLength = 0;
  double gcells = 0;
  double gcells_noseed = 0;
  for (const auto& s : seqs) {
    maxSequenceLength = std::max(maxSequenceLength, s.size());
  }
  for (const auto& c : cmps) {
    gcells += getGcells(c, seedLen, true);
    gcells_noseed += getGcells(c, seedLen, false);
  }

  json stats = {
    {"numComparisons", numComparisons},
    {"maxSequenceLength", maxSequenceLength},
    {"gCells", gcells},
    {"gCellsNoSeed", gcells_noseed},
  };
  return stats;
}

template json getDatasetStats<Comparison>(const ipu::RawSequences& seqs, const std::vector<Comparison>& cmps, int seedLen);
template json getDatasetStats<MultiComparison>(const ipu::RawSequences& seqs, const std::vector<MultiComparison>& cmps, int seedLen);

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

std::vector<std::string> loadStrings(const std::string& path) {
  std::ifstream f(path);
  std::vector<std::string> seqs;
  std::string l;
  while (std::getline(f, l)) seqs.push_back(l);
  return seqs;
}

RawSequences convertStrings(const std::vector<std::string>& seqs) {
  RawSequences ss;
  for (const auto& s : seqs) ss.push_back(s);
  return ss;
}

template<typename C>
SequenceDatabase<C> loadSequences(JsonSequenceConfig& config) {
  PLOGI << "Loading json from " << config.sequences << " " << config.comparisons;
  SequenceDatabase<C> db;
  std::ifstream cmpF(config.comparisons);
  std::ifstream strF(config.sequences);

  #pragma omp parallel sections shared(db)
  {
    #pragma omp section
    db.cmps = json::parse(cmpF);

    #pragma omp section
    db.strings = json::parse(strF).get<std::vector<std::string>>();
  }
  db.seqs = convertStrings(db.strings);
  return std::move(db);
}

template<typename C>
SequenceDatabase<C> loadSequences(SequenceConfig& config, SWConfig& swconfig) {
	SequenceDatabase<C> db;
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
      db.strings.push_back(lH);
      db.strings.push_back(lV);
      const int sizeA = (int) db.strings[(int) (2 * i)].size();
      const int sizeB = (int) db.strings[(int) (2 * i + 1)].size();
      Comparison cmp = {
        i,
        (int) (2 * i),
        sizeA,
        (int) (2 * i + 1),
        sizeB,
        {{{sH1, sV1}, {sH2, sV2}}},
      };
      add_comparison<C>(db.cmps, cmp, swconfig.seedLength);
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
  db.seqs = convertStrings(db.strings);
  return std::move(db);
}

template SequenceDatabase<MultiComparison> loadSequences<MultiComparison>(SequenceConfig& config, SWConfig& swconfig);
template SequenceDatabase<Comparison> loadSequences<Comparison>(SequenceConfig& config, SWConfig& swconfig);

bool isEmpty(GeneratorConfig& c) {
	return c.generatorCount == 0 || c.generatorSeqLen == 0;
}

bool isEmpty(SequenceConfig& c) {
  return c.seqsH.size() == 0 || c.seqsV.size() == 0;
}

bool isEmpty(JsonSequenceConfig& c) {
  return c.sequences.size() == 0 || c.comparisons.size() == 0;
}

SequenceDatabase<Comparison> LoaderConfig::getSequences(ipu::SWConfig swconfig) {
  if (!isEmpty(seqConfig)){
    return loadSequences<Comparison>(seqConfig, swconfig);
  } else if (!isEmpty(jsConfig)) {
    return loadSequences<Comparison>(jsConfig);
  } else if (!isEmpty(genConfig)) {
    return loadSequences<Comparison>(genConfig, swconfig);
  } else {
    throw std::runtime_error("Did not have any sequence info.");
  }
}

SequenceDatabase<MultiComparison> LoaderConfig::getMultiSequences(ipu::SWConfig swconfig) {
  if (!isEmpty(seqConfig)){
    return loadSequences<MultiComparison>(seqConfig, swconfig);
  } else if (!isEmpty(jsConfig)) {
    return loadSequences<MultiComparison>(jsConfig);
  } else if (!isEmpty(genConfig)) {
    return loadSequences<MultiComparison>(genConfig, swconfig);
  } else {
    throw std::runtime_error("Did not have any sequence info.");
  }
}
}