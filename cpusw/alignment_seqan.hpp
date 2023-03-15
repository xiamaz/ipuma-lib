#pragma once

#include <seqan/align.h>
#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

using namespace seqan;

namespace cpu {
  class SeqanAligner {
  public:
    static int align(const std::string_view& seqH, const std::string_view& seqV, int32_t seedHStart, int32_t seedVStart, int32_t seedLen, const ipu::SWConfig& config) {
      Dna5String sH(std::string{seqH});
      Dna5String sV(std::string{seqV});

      Seed<Simple> seed(seedHStart, seedVStart, seedLen);
      Score<int, Simple> scoringScheme(config.matchValue, config.mismatchValue, config.gapExtend);
      auto scoreBoth = extendSeed(seed, sH, sV, seqan::ExtensionDirection::EXTEND_BOTH, scoringScheme, config.xDrop, seedLen, GappedXDrop());

      return scoreBoth;
    }
  };
}