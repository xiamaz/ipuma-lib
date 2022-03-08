#ifndef ALIGNERS_H
#define ALIGNERS_H
#include <thread>
#include <future>

#include <nlohmann/json.hpp>
#include <ssw/ssw_core.hpp>
#include <ssw/ssw.hpp>

#include "cpuswconfig.hpp"
#include "swatlib/timing.hpp"
using json = nlohmann::json;

const static std::vector<int8_t> mat50 = {
 //  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
   	5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,	// A
    -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,	// R
    -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,	// N
    -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,	// D
    -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,	// C
    -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,	// Q
    -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,	// E
   	0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,	// G
    -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,	// H
    -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,	// I
    -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,	// L
    -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,	// K
    -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,	// M
    -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,	// F
    -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,	// P
   	1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,	// S
   	0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, 	// T
    -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, 	// W
    -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, 	// Y
   	0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, 	// V
    -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, 	// B
    -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, 	// Z
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, 	// X
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1 	// *
};

/* This table is used to transform amino acid letters into numbers. */
const static std::vector<int8_t> aa_table = {
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};

struct AlignmentRange {
	int32_t begin;
	int32_t end;
};

struct AlignmentResults {
	std::vector<int32_t> scores;
	std::vector<AlignmentRange> a_range;
	std::vector<AlignmentRange> b_range;

	AlignmentResults(int size) {
		scores.resize(size);
		a_range.resize(size);
		b_range.resize(size);
	}
};

class CPUAligner {
protected:
	CpuSwConfig config;
public:
	CPUAligner(CpuSwConfig c) : config(c) {}

	virtual void compare(const std::vector<std::string>& A, const std::vector<std::string>& B) {
		throw std::runtime_error("Not implemented");
	}
	virtual AlignmentResults get_results() {
		throw std::runtime_error("Not implemented");
	}
};

void testRun() {
	std::cout << "Yell from here\n";
}

class SSWAligner : public CPUAligner {
	std::vector<StripedSmithWaterman::Alignment> alns;
public:
	SSWAligner(CpuSwConfig c) : CPUAligner(c) {}

	void ssw_aligner_align(const StripedSmithWaterman::Aligner& ssw_aligner, const char* query, const int query_len, const char* ref, const int ref_len, const StripedSmithWaterman::Filter& filter, StripedSmithWaterman::Alignment& alignment, const int32_t maskLen, swatlib::TickTock& t) {
    int8_t* translated_query = new int8_t[query_len];
    ssw_aligner.TranslateBase(query, query_len, translated_query);

    // calculate the valid length
    int valid_ref_len = ref_len;
    int8_t* translated_ref = new int8_t[valid_ref_len];
    ssw_aligner.TranslateBase(ref, valid_ref_len, translated_ref);


    const int8_t score_size = 2;
    s_profile* profile = ssw_init(translated_query, query_len, ssw_aligner.score_matrix_,
                                  ssw_aligner.score_matrix_size_, score_size);

    uint8_t flag = 0;
    if (filter.report_begin_position) flag |= 0x08;
    if (filter.report_cigar) flag |= 0x0f;

		t.tick();
    s_align* s_al = ssw_align(profile, translated_ref, valid_ref_len,
                              static_cast<int>(ssw_aligner.gap_opening_penalty_),
                              static_cast<int>(ssw_aligner.gap_extending_penalty_),
                              flag, filter.score_filter, filter.distance_filter, maskLen);
		t.tock();

    alignment.Clear();
    alignment.sw_score           = s_al->score1;
    alignment.ref_begin          = s_al->ref_begin1;
    alignment.ref_end            = s_al->ref_end1;
    alignment.query_begin        = s_al->read_begin1;
    alignment.query_end          = s_al->read_end1;
    // alignment->mismatches = CalculateNumberMismatch(alignment, translated_ref, translated_query, query_len);

    // Free memory
    delete [] translated_query;
    delete [] translated_ref;
    align_destroy(s_al);
    init_destroy(profile);
	}

	void compare_thread(int workerId, const std::vector<std::string>& A, const std::vector<std::string>& B, swatlib::TickTock& t) {
		// calculate comparison
    // int matSize = static_cast<int>(std::sqrt(mat50.size()));
		StripedSmithWaterman::Aligner ssw_aligner; //(mat50.data(), matSize, aa_table.data(), aa_table.size());
		
		if (config.swconfig.datatype == swatlib::DataType::aminoAcid) {
    	int matSize = static_cast<int>(std::sqrt(mat50.size()));
    	ssw_aligner.ReBuild(mat50.data(), matSize, aa_table.data(), aa_table.size());
		} else {
			ssw_aligner.ReBuild(
				config.swconfig.matchValue,
				-config.swconfig.mismatchValue,
				(-config.swconfig.gapInit)-config.swconfig.gapExtend,
				-config.swconfig.gapExtend,
				-config.swconfig.ambiguityValue
			);
		}
		StripedSmithWaterman::Filter ssw_filter(false, false, 0, 32767);
		for (int c = workerId; c < A.size(); c += config.algoconfig.threads) {
			const auto& a = A[c];
			const auto& b = B[c];
			// t.tick();
			ssw_aligner_align(ssw_aligner, a.data(), a.size(), b.data(), b.size(), ssw_filter, alns[c], std::max((int)(b.size() / 2), 15), t);
			// t.tock();
		}
	}

	void compare(const std::vector<std::string>& A, const std::vector<std::string>& B) override {
		swatlib::TickTock outer;
		std::vector<swatlib::TickTock> innerTimers(config.algoconfig.threads);
		std::deque<std::thread> threads;
		alns.resize(A.size());
		outer.tick();
		for (int n = 0; n < config.algoconfig.threads; ++n) {
			threads.push_back(std::thread(&SSWAligner::compare_thread, this, n, std::cref(A), std::cref(B), std::ref(innerTimers[n])));
		}

		for (int n = 0; n < config.algoconfig.threads; ++n) {
			threads[n].join();
		}
		outer.tock();

		double cellCount = 0;
		for (int i = 0; i < A.size(); ++i) {
			cellCount += static_cast<double>(A[i].size() * B[i].size());
		}

		double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
		double innerTime = 0;
		for (const auto& t : innerTimers) {
			innerTime = std::max(innerTime, static_cast<double>(t.accumulate_microseconds()) / 1e6);
		}
		double gcupsOuter = cellCount / outerTime / 1e9;
		double gcupsInner = cellCount / innerTime / 1e9;
		json jsonlog = {
			{"time_outer_s", outerTime},
			{"time_inner_s", innerTime},
			{"gcups_outer", gcupsOuter},
			{"gcups_inner", gcupsInner},
			{"cells", cellCount},
		};
		PLOGW << "CPUJSONLOG" << jsonlog.dump();
	}

	AlignmentResults get_results() override {
		AlignmentResults results(alns.size());
		for (int i = 0; i < alns.size(); ++i) {
			const auto& aln = alns[i];
			results.scores[i] = aln.sw_score;
			results.a_range[i].begin = aln.ref_begin;
			results.a_range[i].end = aln.ref_end;
			results.b_range[i].begin = aln.query_begin;
			results.b_range[i].end = aln.query_end;
		}
		return results;
	}
};

#endif