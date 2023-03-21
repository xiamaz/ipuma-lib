#pragma once

#include <string_view>
#include "ipu_config.h"

extern "C" {
  #include "gaba.h"										/* just include gaba.h */
}


namespace cpu {

  std::tuple<std::vector<uint8_t>, std::vector<uint8_t>> createSeedParts(const std::string_view& sequence, int seedStart, int seedLength, swatlib::Encoding& enc) {
    auto seqRight = sequence.substr(seedStart + seedLength, sequence.size() - (seedStart + seedLength));
    auto seqLeft = sequence.substr(0, seedStart);

    /*
    if (seedStart == 0 || sequence.size() - (seedStart + seedLength) < 1) {
      PLOGF << "seedStart " << seedStart << " size " << sequence.size() << " seedLen " << seedLength;
      exit(1);
    }
    */

    std::string strRight(seqRight);
    std::string strLeft(seqLeft.rbegin(), seqLeft.rend());

    return {
      enc.encode(strLeft),
      enc.encode(strRight),
    };
  }

  int alignGabaInner( const std::vector<uint8_t>& encSeqH, const std::vector<uint8_t>& encSeqV, const gaba_params_t& params, gaba_t* ctx) {
	    char const t[64] = { 0 };							/* tail array */

	    struct gaba_section_s asec = gaba_build_section(0, encSeqH.data(), (uint32_t) encSeqH.size());
	    struct gaba_section_s bsec = gaba_build_section(2, encSeqV.data(), (uint32_t) encSeqV.size());
	    struct gaba_section_s tail = gaba_build_section(4, t, 64);

	    /* create thread-local object */
	    gaba_dp_t *dp = gaba_dp_init(ctx);					/* dp[0] holds a 64-cell-wide context */
	    // gaba_dp_t *dp_32 = &dp[_dp_ctx_index(32)];			/* dp[1] and dp[2] are narrower ones */
	    // gaba_dp_t *dp_16 = &dp[_dp_ctx_index(16)];

	    /* init section pointers */
	    struct gaba_section_s const *ap = &asec, *bp = &bsec;
	    struct gaba_fill_s const *f = gaba_dp_fill_root(dp,	/* dp -> &dp[_dp_ctx_index(band_width)] makes the band width selectable */
	    	ap, 0,											/* a-side (reference side) sequence and start position */
	    	bp, 0,											/* b-side (query) */
	    	UINT32_MAX										/* max extension length */
	    );

	    /* until X-drop condition is detected */
	    struct gaba_fill_s const *m = f;					/* track max */
	    while((f->status & GABA_TERM) == 0) {
	    	if(f->status & GABA_UPDATE_A) { ap = &tail; }	/* substitute the pointer by the tail section's if it reached the end */
	    	if(f->status & GABA_UPDATE_B) { bp = &tail; }

	    	f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);	/* extend the banded matrix */
	    	m = f->max > m->max ? f : m;					/* swap if maximum score was updated */
	    }

      auto score = m->max;

	    /* clean up */
	    gaba_dp_clean(dp);
      return score;
  }

  class GabaAligner {
    gaba_t* ctx = nullptr;
    gaba_params_t params;

    const ipu::SWConfig& config;
  public:
    GabaAligner(const ipu::SWConfig& config) : config(config) {
      int8_t _m = (int8_t) config.matchValue;
      int8_t _x = (int8_t) config.mismatchValue;
      params = {
        .score_matrix = {
          _m, _x, _x, _x,
          _x, _m, _x, _x,
          _x, _x, _m, _x,
          _x, _x, _x, _m,
        },
        .gi = (int8_t) -config.gapInit,
        .ge = (int8_t) -config.gapExtend,
        .xdrop = (int8_t) config.xDrop,
      };
      ctx = gaba_init(&params);
    }

    ~GabaAligner() {
      gaba_clean(ctx);
    }

    int align(const std::string_view& seqH, const std::string_view& seqV, int32_t seedHStart, int32_t seedVStart, int32_t seedLen) {
      swatlib::Encoding enc({
        {'A', 0x01},
        {'C', 0x02},
        {'G', 0x04},
        {'T', 0x08},
      });
      auto [seqHL, seqHR] = createSeedParts(seqH, seedHStart, seedLen, enc);
      auto [seqVL, seqVR] = createSeedParts(seqV, seedVStart, seedLen, enc);

      int scoreLeft = 0, scoreRight = 0;
      if (seedHStart > 0 && seedVStart > 0) {
        scoreLeft = alignGabaInner(seqHL, seqVL, params, ctx);
      }
      if (seedHStart + seedLen < seqH.size() && seedVStart + seedLen < seqV.size()) {
        scoreRight = alignGabaInner(seqHR, seqVR, params, ctx);
      }

      return scoreLeft + scoreRight + seedLen;
    }
  };
}