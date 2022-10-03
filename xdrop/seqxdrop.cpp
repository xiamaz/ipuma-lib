// Source: https://docs.seqan.de/seqan/3-master-user/tutorial_pairwise_alignment.html
#include <utility>
 
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "./data.h"

seqan3::dna4_vector convertSequence(const std::string& string) {
    seqan3::dna4_vector dna4_str{};
    for (auto c : string) dna4_str.push_back(seqan3::assign_char_to(c, seqan3::dna4{}));
    return dna4_str;
}

std::vector<seqan3::dna4_vector> convertSequences(const std::vector<std::string>& strings) {
    std::vector<seqan3::dna4_vector> seqs{};
    for (const auto& s : strings) {
        seqs.push_back(convertSequence(s));
    }
    return seqs;
}

std::vector<int> seqanAlign(const std::vector<std::string>& queryStrs, const std::vector<std::string>& referenceStrs) {
    auto queries = convertSequences(queryStrs);
    auto references = convertSequences(referenceStrs);
    // seqan3::dna4_vector query = convertSequence(queryStr);
    // seqan3::dna4_vector reference = convertSequence(referenceStr);

    // Configure the alignment kernel.
    auto config =
        seqan3::align_cfg::method_global(
            seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
            seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}
        ) | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-1}};;

    auto results = seqan3::align_pairwise(seqan3::views::zip(queries, references), config);
    std::vector<int> scores{};
    for (auto const& res : results) {
        scores.push_back(res.score());
    }
    return scores;
}

int main()
{
    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    auto score = seqanAlign(TEST_refs, TEST_queries);
    seqan3::debug_stream << "Score: " << score << '\n';
}