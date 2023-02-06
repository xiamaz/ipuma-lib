#pragma once

#include <cmath>
#include "gtest/gtest.h"

#include "ssw/ssw.hpp"

#include <seqan/align.h>
#include <seqan/stream.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

int alignFullSSW(const std::string cpp_seqH, const std::string cpp_seqV) {
    StripedSmithWaterman::Aligner ssw_aligner(1, 1, 1, 1, 1);

    StripedSmithWaterman::Filter ssw_filter;
    StripedSmithWaterman::Alignment aln;
    ssw_aligner.Align(cpp_seqH.data(), cpp_seqH.size(), cpp_seqV.data(), cpp_seqV.size(), ssw_filter, &aln, std::max((int)(cpp_seqH.size() / 2), 15));
    return aln.sw_score;
}

int alignSeqSeqan(const std::string cpp_seqH, const std::string cpp_seqV, const int XDrop) {
    using namespace seqan;
    // The horizontal and vertical sequence (subject and query sequences).
    Dna5String seqH = cpp_seqH;
    Dna5String seqV = cpp_seqV;
    // Create the seed sequence.
    Seed<Simple> seed(0, 0, 0, 0);

    // Perform match extension.
    Score<int, Simple> scoringScheme(1, -1, -1);
    auto score = extendSeed(seed, seqH, seqV, EXTEND_RIGHT, scoringScheme, XDrop, GappedXDrop());

    if (true) {
      // Perform a banded alignment.
      Align<Dna5String> align;
      resize(rows(align), 2);
      assignSource(row(align, 0), infix(seqH, beginPositionH(seed),
                                        endPositionH(seed)));
      assignSource(row(align, 1), infix(seqV, beginPositionV(seed),
                                        endPositionV(seed)));

      globalAlignment(align, scoringScheme);
    
      std::cout << "Resulting alignment\n" << align << "\n";
      std::cout << "\n";
    }
    return score;
}

inline std::string aln2string(StripedSmithWaterman::Alignment &aln) {
  std::stringstream ss;
  ss << "score=" << aln.sw_score << " score2=" << aln.sw_score_next_best;
  ss << " rbegin=" << aln.ref_begin << " rend=" << aln.ref_end;
  ss << " qbegin=" << aln.query_begin << " qend=" << aln.query_end;
  ss << " rend2=" << aln.ref_end_next_best << " mismatches=" << aln.mismatches;
  ss << " cigarstr=" << aln.cigar_string;
  return ss.str();
}

std::string initerstring(int minlen) {
  std::string x = "";
  srand(0);
  std::vector bp{"C", "A", "T", "G"};
  for (int i = 0; i < minlen; i++) {
    x += bp[rand()%4];
  }
  return x;
}