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
      auto scoreRight = extendSeed(seed, sH, sV, EXTEND_RIGHT, scoringScheme, config.xDrop, GappedXDrop());
      auto scoreLeft = extendSeed(seed, sH, sV, EXTEND_LEFT, scoringScheme, config.xDrop, GappedXDrop());
      PLOGE << scoreRight;
      PLOGE << scoreLeft;
      return scoreRight + scoreLeft;
    }
  };
}