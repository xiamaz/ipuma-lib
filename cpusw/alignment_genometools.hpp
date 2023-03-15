#pragma once

#include "genometools.h"
#include "match/diagbandseed.h"
#include "match/xdrop.h"

#include "ipu_config.h"

namespace cpu {

  class GenomeToolsAligner {
  public:
    static int align(const std::string_view& seqH, const std::string_view& seqV, int32_t seedHStart, int32_t seedVStart, int32_t seedLen, const ipu::SWConfig& config) {
      gt_lib_init();
      GtAlphabet* alpha = gt_alphabet_new_dna();

      GtEncseqBuilder* eb = gt_encseq_builder_new(alpha);
      gt_encseq_builder_disable_description_support(eb);

      GtSeqabstract *useq, *vseq;
      useq = gt_seqabstract_new_gtuchar(
        true, GT_READMODE_FORWARD, (const GtUchar*) seqH.data(), seqH.size(), 0, seqH.size());
      vseq = gt_seqabstract_new_gtuchar(
        true, GT_READMODE_FORWARD, (const GtUchar*) seqV.data(), seqV.size(), 0, seqV.size());

      GtXdropArbitraryscores xdropConfig{
        .mat = config.matchValue,
        .mis = config.mismatchValue,
        .ins = config.gapExtend,
        .del = config.gapExtend,
      };

      auto* xres = gt_xdrop_resources_new(&xdropConfig);

      GtXdropbest result;

      PLOGF << "GTOOLS " << result.score;

      int score;

      gt_evalxdroparbitscoresextend(
        true, // forward
        &result,
        xres,
        useq,
        vseq,
        config.xDrop
      );

      score = result.score;

      gt_seqabstract_delete(useq);
      gt_seqabstract_delete(vseq);

      gt_xdrop_resources_delete(xres);
      gt_encseq_builder_delete(eb);
      gt_alphabet_delete(alpha);

      gt_lib_clean();

      return score;
    }
  };
}