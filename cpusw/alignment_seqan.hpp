#pragma once

#include <seqan/align.h>
#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

using namespace seqan;

void align(const std::string& seqH, const std::string& seqV, int32_t seedHStart, int32_t seedVStart) {
  Dna5String sH = seqH;
  Dna5String sV = seqV;

  Seed<Simple> seed(seedHStart, seedVStart, 17);
  Score<int, Simple> scoringScheme(1, -1, -1);
  auto scoreRight = extendSeed(seed, sH, sV, EXTEND_RIGHT, scoringScheme, 10, GappedXDrop());
  auto scoreLeft = extendSeed(seed, sH, sV, EXTEND_LEFT, scoringScheme, 10, GappedXDrop());
}

void runAlignment(const ipu::RawSequences& seqs, const ipu::Comparisons& cmps) {
  swatlib::TickTock t;
  t.tick();
  double cells = 0;
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
    align(seqs[cmp.indexA], seqs[cmp.indexB], cmp.seedAStartPos, cmp.seedBStartPos);
  }
  t.tock();
  double gcups = cells / t.accumulate_microseconds() * 1e6;
  PLOGE << gcups;
}