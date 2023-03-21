#pragma once

#include <string_view>
#include "ipu_config.h"
#include "ksw2.h"

namespace cpu {
  int alignKsw2(const std::vector<uint8_t>& H, const std::vector<uint8_t>& V, const ipu::SWConfig& config) {
      ksw_extz_t result;
      int bandwidth = -1;
      int zdrop = config.xDrop; // default should be Z-Drop = 400
      int8_t a = config.matchValue;
      int8_t b = config.mismatchValue;
      int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };

      result.max_q = result.max_t = result.mqe_t = result.mte_q = -1;
	    result.max = 0, result.mqe = result.mte = KSW_NEG_INF;
	    result.n_cigar = 0;

      ksw_extz2_sse(
        nullptr,
        H.size(),
        H.data(),
        V.size(),
        V.data(),
        5,
        mat,
        -config.gapInit,
        -config.gapExtend,
        bandwidth,
        zdrop,
        0,
        // KSW_EZ_SCORE_ONLY,
        KSW_EZ_SCORE_ONLY | KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP | KSW_EZ_EXTZ_ONLY,
        &result
      );

      return result.max;
  }

  class Ksw2Aligner {
    const ipu::SWConfig& config;
    swatlib::Encoding enc;
  public:
    Ksw2Aligner(const ipu::SWConfig& config) : config(config), enc({
        {'A', 0x00},
        {'C', 0x01},
        {'G', 0x02},
        {'T', 0x03},
    }) {}

    int align(const std::string_view& seqH, const std::string_view& seqV, int32_t seedHStart, int32_t seedVStart, int32_t seedLen) {
      auto [seqHL, seqHR] = createSeedParts(seqH, seedHStart, seedLen, enc);
      auto [seqVL, seqVR] = createSeedParts(seqV, seedVStart, seedLen, enc);

      int scoreLeft = 0, scoreRight = 0;
      if (seedHStart > 0 && seedVStart > 0) {
        scoreLeft = alignKsw2(seqHL, seqVL, config);
      }
      if (seedHStart + seedLen < seqH.size() && seedVStart + seedLen < seqV.size()) {
        scoreRight = alignKsw2(seqHR, seqVR, config);
      }

      return scoreLeft + scoreRight + seedLen;
    }
  };
}