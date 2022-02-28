#ifndef ALIGNERS_H
#define ALIGNERS_H
#include <thread>
#include <future>

#include <ssw/ssw_core.hpp>
#include <ssw/ssw.hpp>

#include "cpuswconfig.hpp"
#include "swatlib/timing.hpp"

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
		StripedSmithWaterman::Aligner ssw_aligner(
			config.swconfig.matchValue,
			-config.swconfig.mismatchValue,
			(-config.swconfig.gapInit)-config.swconfig.gapExtend,
			-config.swconfig.gapExtend,
			-config.swconfig.ambiguityValue
		);
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
		PLOGI << "Cell Count: " << cellCount;

		double outerTime = static_cast<double>(outer.accumulate_microseconds()) / 1e6;
		double innerTime = 0;
		for (const auto& t : innerTimers) {
			innerTime = std::max(innerTime, static_cast<double>(t.accumulate_microseconds()) / 1e6);
		}
		PLOGI << "outer: " << outerTime << "s inner: " << innerTime << "s";
		double gcupsOuter = cellCount / outerTime / 1e9;
		double gcupsInner = cellCount / innerTime / 1e9;
		PLOGI << "GCUPS outer: " << gcupsOuter << " inner: " << gcupsInner;
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