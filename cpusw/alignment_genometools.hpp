#pragma once

#include "genometools.h"
#include "match/diagbandseed.h"

void test() {
  gt_lib_init();
  GtAlphabet* alpha = gt_alphabet_new_dna();

  GtEncseqBuilder* eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_disable_description_support(eb);

  std::string teststr1 = "aaaagggttttt";
  std::string teststr2 = "caaagggtttgt";

  gt_encseq_builder_add_cstr(eb, teststr1.c_str(), teststr2.size(), NULL);
  gt_encseq_builder_add_cstr(eb, teststr2.c_str(), teststr2.size(), NULL);

  GtEncseq* es = gt_encseq_builder_build(eb, NULL);
  
GtDiagbandseedInfo *gt_diagbandseed_info_new(es,
                                             es,
                                             5,
                                             42,
                                             unsigned int spacedseedweight,
                                             unsigned int seedlength,
                                             bool norev,
                                             bool nofwd,
                                             const GtRange *seedpairdistance,
                                             GtDiagbandseedBaseListType splt,
                                             GtDiagbandseedBaseListType kmplt,
                                             bool verify,
                                             bool verbose,
                                             bool debug_kmer,
                                             bool debug_seedpair,
                                             bool use_kmerfile,
                                             bool trimstat_on,
                                             GtUword maxmat,
                                             const GtStr *chainarguments,
                                             const GtStr
                                               *diagband_statistics_arg,
                                             const GtDiagbandseedExtendParams
                                               *extp);

  gt_encseq_delete(es);
  gt_encseq_builder_delete(eb);
  gt_alphabet_delete(alpha);

  gt_lib_clean();
}