// Source: https://docs.seqan.de/seqan/3-master-user/tutorial_pairwise_alignment.html
#include <stdlib.h>
#include <cstring>
#include <utility>
 
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "../src/swatlib/swatlib.h"
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
        ) | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{2}, seqan3::mismatch_score{-2}}}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-3}};;

    auto results = seqan3::align_pairwise(seqan3::views::zip(queries, references), config);
    std::vector<int> scores{};
    for (auto const& res : results) {
        scores.push_back(res.score());
    }
    return scores;
}

int xdropAlign(const std::string& s1, const std::string& s2) {
    int N = s2.size();
    int M = s1.size();

    swatlib::Matrix<int> H(M + 1, N + 1, 0); // DEBUG

    int mis = -2;
    int mat = 2;
    int X = 10;

    int inf = 10000;

    auto calcSPrime = [&](int ij, int dis){
        return std::floor(
            ij * ((double) mat /2) - dis * (mat - mis)
        );
    };

    int i = 0;
    for (; i < std::min(M, N) && s1[i] == s2[i]; ++i){
            H(i, i) = H(i, i)+1; // DEBUG
    }

    int R_offset = N;

    int* R0 = (int*) malloc(2 * N * sizeof(int));
    int* R1 = (int*) malloc(2 * N * sizeof(int));
    memset(R0, -inf, 2*N*sizeof(int));
    memset(R1, -inf, 2*N*sizeof(int));
    R0[R_offset] = i;

    int xdrop_offset = ((X + mat / 2) / (mat - mis)) + 1;
    int* T = (int*) malloc((M+N+xdrop_offset+1)*sizeof(int));
    memset(T, 0, (M+N+xdrop_offset+1)*sizeof(int));
    T[xdrop_offset] = calcSPrime(i + i, 0);
    int Tt = T[xdrop_offset];

    int d = 0;
    int Lmin = 0;
    int Lmax = -inf;
    int L = 0;
    int Umin = inf;
    int Umax = 0;
    int U = 0;

    if (i == N) Lmax = 0;
    if (i == M) Umin = 0;

    for (int d = 1; d < M + N; ++d) {
        int dd = d - xdrop_offset;
        for (int k = L-1; k <= U+1; ++k) {
            int i = -1;
            if (L < k) {
                i = std::max(i, R0[k-1+R_offset]+1);
            }
            if (L <= k && k <= U) {
                i = std::max(i, R0[k+R_offset]+1);
            }
            if (k < U) {
                i = std::max(i, R0[k+1+R_offset]);
            }
            int j = i - k;
            int cs = calcSPrime(i+j, d);
            if (i >= 0 && j>=0 && i < M && j < N) H(i, j) = H(i, j)+1; // DEBUG
            // std::cout << i << " " << j << " " << cs << std::endl;
            if (i >= 0 && j >= 0 && cs >= T[dd+xdrop_offset] - X) {
                while (i < M && j < N && s1[i] == s2[j]) {
                    H(i, j) = H(i, j)+1; // DEBUG
                    ++i, ++j;
                }
                R1[k + R_offset] = i;
                std::cout << i << " ,j=" << j << " ,M=" << M << " ,N=" << N << " ,cs=" << cs << std::endl;
                cs = calcSPrime(i+j, d);
                Tt = std::max(Tt, cs);
            } else {
                R1[k + R_offset] = -inf;
            }
        }
        // std::cout << Tt << std::endl;
        T[d+xdrop_offset] = Tt;

        for (int rk = 0; rk < 2*N; ++rk) {
            int i = R1[rk];
            int k = rk - R_offset;
            if (i == N+k) Lmax = std::max(Lmax, k);
            if (i == M) Lmax = std::min(Umin, k);
            if (i > -inf) {
                Umax = std::max(Umax, k);
                Lmin = std::min(Lmin, k);
            }
        }
        L = std::max(Lmin, Lmax+2);
        U = std::min(Umax, Umin-2);

        int* h = R1;
        R1 = R0;
        R0 = h;

        if (L > U + 2) break;
    }

    free(R0);
    free(R1);
    free(T);
    std::cout << H.toString();
    return Tt;
}

int main()
{
    for (int i = 0; i < TEST_refs.size(); ++i) {
        auto score = xdropAlign(TEST_refs[i], TEST_queries[i]);
        seqan3::debug_stream << "s1=" << TEST_refs[i] << "\ns2=" << TEST_queries[i] << "\n" << "   Score=" << score << "\n######\n";
    }
    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    // auto score = seqanAlign(TEST_refs, TEST_queries);
    // for (int i = 0; i < TEST_refs.size(); ++i) {
    //     seqan3::debug_stream << TEST_refs[i] << " " << TEST_queries[i] << " " << "Score: " << score[i] << '\n';
    // }
}