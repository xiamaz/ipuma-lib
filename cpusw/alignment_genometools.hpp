#pragma once

#include "genometools.h"
#include "match/diagbandseed.h"
#include "match/xdrop.h"

#include "ipu_config.h"

namespace cpu {

  class GenomeToolsAligner {
  public:
    static int align(const std::string_view& seqH, const std::string_view& seqV, int32_t seedHStart, int32_t seedVStart, int32_t seedLen, const ipu::SWConfig& config) {

      GtSeqabstract *useq_right, *vseq_right;
      GtSeqabstract *useq_left, *vseq_left;
      useq_left = gt_seqabstract_new_gtuchar(
        false, GT_READMODE_FORWARD, (const GtUchar*) seqH.data(), seedHStart, 0, seqH.size());
      vseq_left = gt_seqabstract_new_gtuchar(
        false, GT_READMODE_FORWARD, (const GtUchar*) seqV.data(), seedVStart, 0, seqV.size());
      useq_right = gt_seqabstract_new_gtuchar(
        true, GT_READMODE_FORWARD, (const GtUchar*) seqH.data(), seqH.size() - (seedHStart + seedLen), seedHStart + seedLen, seqH.size());
      vseq_right = gt_seqabstract_new_gtuchar(
        true, GT_READMODE_FORWARD, (const GtUchar*) seqV.data(), seqV.size() - (seedVStart + seedLen), seedVStart + seedLen, seqV.size());

      GtXdropArbitraryscores xdropConfig{
        // .mat = config.matchValue,
        // .mis = config.mismatchValue,
        // .ins = config.gapExtend,
        // .del = config.gapExtend,
        .mat = 2,
        .mis = -1,
        .ins = -2,
        .del = -2,
      };

      auto* xres_left = gt_xdrop_resources_new(&xdropConfig);
      auto* xres_right = gt_xdrop_resources_new(&xdropConfig);

      GtXdropbest result_left, result_right;

      int score;

      gt_evalxdroparbitscoresextend(
        true, // forward
        &result_left,
        xres_left,
        useq_left,
        vseq_left,
        config.xDrop
      );

      gt_evalxdroparbitscoresextend(
        true, // forward
        &result_right,
        xres_right,
        useq_right,
        vseq_right,
        config.xDrop
      );

      score = result_left.score + result_right.score;

      gt_seqabstract_delete(useq_left);
      gt_seqabstract_delete(vseq_left);
      gt_seqabstract_delete(useq_right);
      gt_seqabstract_delete(vseq_right);

      gt_xdrop_resources_delete(xres_left);
      gt_xdrop_resources_delete(xres_right);

      return score;
    }
  };
}