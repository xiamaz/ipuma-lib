// ssw_cpp.cpp
// Created by Wan-Ping Lee
// Last revision by Mengyao Zhao on 2017-05-30

#include "ssw.hpp"
#include "ssw_core.hpp"

#include <sstream>
#include <cassert>

namespace {

  static const int8_t kBaseTranslation[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    //   A     C            G
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    //             T
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    //   a     c            g
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    //             t
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
  };

  void BuildSwScoreMatrix(const uint8_t& match_score,
                          const uint8_t& mismatch_penalty,
                          const uint8_t& ambiguity_penalty,
                          int8_t* matrix) {

    // The score matrix looks like
    //                 // A,  C,  G,  T,  N
    //  score_matrix_ = { 2, -2, -2, -2, -2, // A
    //                   -2,  2, -2, -2, -2, // C
    //                   -2, -2,  2, -2, -2, // G
    //                   -2, -2, -2,  2, -2, // T
    //                   -2, -2, -2, -2, -2};// N

    int id = 0;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        matrix[id] = ((i == j) ? match_score : static_cast<int8_t>(-mismatch_penalty));
        ++id;
      }
//    matrix[id] = static_cast<int8_t>(-mismatch_penalty); // For N
      matrix[id] = static_cast<int8_t>(-ambiguity_penalty); // For N
      ++id;
    }

    for (int i = 0; i < 5; ++i)
//    matrix[id++] = static_cast<int8_t>(-mismatch_penalty); // For N
      matrix[id++] = static_cast<int8_t>(-ambiguity_penalty); // For N
  }

  void ConvertAlignment(const s_align& s_al,
                        const int& query_len,
                        StripedSmithWaterman::Alignment* al) {
    al->sw_score           = s_al.score1;
    al->sw_score_next_best = s_al.score2;
    al->ref_begin          = s_al.ref_begin1;
    al->ref_end            = s_al.ref_end1;
    al->query_begin        = s_al.read_begin1;
    al->query_end          = s_al.read_end1;
    al->ref_end_next_best  = s_al.ref_end2;

    al->cigar.clear();
    al->cigar_string.clear();

    if (s_al.cigarLen > 0) {
      std::ostringstream cigar_string;
      if (al->query_begin > 0) {
        uint32_t cigar = to_cigar_int(al->query_begin, 'S');
        al->cigar.push_back(cigar);
        cigar_string << al->query_begin << 'S';
      }

      for (int i = 0; i < s_al.cigarLen; ++i) {
        al->cigar.push_back(s_al.cigar[i]);
        cigar_string << cigar_int_to_len(s_al.cigar[i]) << cigar_int_to_op(s_al.cigar[i]);
      }

      int end = query_len - al->query_end - 1;
      if (end > 0) {
        uint32_t cigar = to_cigar_int(end, 'S');
        al->cigar.push_back(cigar);
        cigar_string << end << 'S';
      }

      al->cigar_string = cigar_string.str();
    } // end if
  }

// @Function:
//     Calculate the length of the previous cigar operator
//     and store it in new_cigar and new_cigar_string.
//     Clean up in_M (false), in_X (false), length_M (0), and length_X(0).
  void CleanPreviousMOperator(
    bool* in_M,
    bool* in_X,
    uint32_t* length_M,
    uint32_t* length_X,
    std::vector<uint32_t>* new_cigar,
    std::ostringstream* new_cigar_string) {
    if (*in_M) {
      uint32_t match = to_cigar_int(*length_M, '=');
      new_cigar->push_back(match);
      (*new_cigar_string) << *length_M << '=';
    } else if (*in_X) { //in_X
      uint32_t match = to_cigar_int(*length_X, 'X');
      new_cigar->push_back(match);
      (*new_cigar_string) << *length_X << 'X';
    }

    // Clean up
    *in_M = false;
    *in_X = false;
    *length_M = 0;
    *length_X = 0;
  }

// @Function:
//     1. Calculate the number of mismatches.
//     2. Modify the cigar string: (if NO_SAM_REPLACE_M_WITH_MATCH_MISMATCH is not set)
//         differentiate alignment matches (M) into sequence match (=) and sequence mismatch (X).
// @Return:
//     The number of mismatches.
  int CalculateNumberMismatch(
    StripedSmithWaterman::Alignment* al,
    int8_t const *ref,
    int8_t const *query,
    const int& query_len) {

    if (al->cigar.empty()) return 0;

    ref   += al->ref_begin;
    
    int mismatch_length = 0;

    std::vector<uint32_t> new_cigar;
    std::ostringstream new_cigar_string;
    
    int start_op = 0;
    if (al->query_begin > 0) {
      uint32_t cint = to_cigar_int(al->query_begin, 'S');
      if (cint != al->cigar[0]) {
        char op = cigar_int_to_op(al->cigar[0]);
        if (op == 'S') {
          assert(false && "This should not happen! First CIGAR operation is soft clip but is the wrong value!");
          abort();
        }
        // needs a new Soft clip entry
        new_cigar.push_back(cint);
        new_cigar_string << al->query_begin << 'S';
        query += al->query_begin;
      }
    } // else it will be processed normally
    

    bool in_M = false; // the previous is match
    bool in_X = false; // the previous is mismatch
    uint32_t length_M = 0;
    uint32_t length_X = 0;

    for (unsigned int i = start_op; i < al->cigar.size(); ++i) {
      char op = cigar_int_to_op(al->cigar[i]);
      uint32_t length = cigar_int_to_len(al->cigar[i]);
      if (op == 'M') {
#ifndef NO_SAM_REPLACE_M_WITH_MATCH_MISMATCH
        // replace the old M op
        for (uint32_t j = 0; j < length; ++j) {
          if (*ref != *query) {
            ++mismatch_length;
            if (in_M) { // the previous is match; however the current one is mismatch
              CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
            }
            // start X op
            length_M = 0;
            ++length_X;
            in_M = false;
            in_X = true;
          } else { // *ref == *query
            if (in_X) { // the previous is mismatch; however the current one is match
              CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
            }
            // start = op
            ++length_M;
            length_X = 0;
            in_M = true;
            in_X = false;
          } // end of if (*ref != *query)
          ++ref;
          ++query;
        }
#else
        assert(!in_X && "No mismatch possible in 'M' cigar op");
        // keep the old M op
        query += length;
        ref += length;
        in_M = true;
        length_M += length;
#endif
      } else if (op == 'I') {
        // add an I op
        query += length;
        mismatch_length += length;
        CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
        new_cigar.push_back(al->cigar[i]);
        new_cigar_string << length << op;
      } else if (op == 'D') {
        // add a D op
        ref += length;
        mismatch_length += length;
        CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
        new_cigar.push_back(al->cigar[i]);
        new_cigar_string << length << op;
      } else {
        assert((op == 'X' || op == '=' || op == 'S') && "No unexpected op"); // no support for N, H, or P yet
        if (op == 'S') {
          query += length;
        } else if (op == 'X' || op == '=') {
          // X or =
          query += length;
          ref += length;
        } else {
          abort();
        }
        CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
        new_cigar.push_back(al->cigar[i]);
        new_cigar_string << length << op;
      }
    }

    CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);

    int end = query_len - al->query_end - 1;
    if (end > 0) {
      uint32_t cint = to_cigar_int(end, 'S');
      if (cint != al->cigar[al->cigar.size()-1]) {
        char op = cigar_int_to_op(al->cigar[al->cigar.size()-1]);
        if (op == 'S') {
          assert(false && "This should not happen! Last CIGAR operation is soft clip but is the wrong value!");
          abort();
        }
        new_cigar.push_back(cint);
        new_cigar_string << end << 'S';
      }
    }

    al->cigar_string = new_cigar_string.str();
    al->cigar.swap(new_cigar);
    return mismatch_length;
  }

  void SetFlag(const StripedSmithWaterman::Filter& filter, uint8_t* flag) {
    if (filter.report_begin_position) *flag |= 0x08;
    if (filter.report_cigar) *flag |= 0x0f;
  }

// http://www.cplusplus.com/faq/sequences/arrays/sizeof-array/#cpp
  template <typename T, size_t N>
  inline size_t SizeOfArray( const T(&)[ N ] )
  {
    return N;
  }

} // namespace



namespace StripedSmithWaterman {

  Aligner::Aligner(void)
    : score_matrix_(NULL)
    , score_matrix_size_(5)
    , translation_matrix_(NULL)
    , match_score_(2)
    , mismatch_penalty_(2)
    , gap_opening_penalty_(3)
    , gap_extending_penalty_(1)
    , ambiguity_penalty_(2)
    , translated_reference_(NULL)
    , reference_length_(0)
  {
    BuildDefaultMatrix();
  }

  Aligner::Aligner(
    const uint8_t& match_score,
    const uint8_t& mismatch_penalty,
    const uint8_t& gap_opening_penalty,
    const uint8_t& gap_extending_penalty,
    const uint8_t& ambiguity_penalty)

    : score_matrix_(NULL)
    , score_matrix_size_(5)
    , translation_matrix_(NULL)
    , match_score_(match_score)
    , mismatch_penalty_(mismatch_penalty)
    , gap_opening_penalty_(gap_opening_penalty)
    , gap_extending_penalty_(gap_extending_penalty)
    , ambiguity_penalty_(ambiguity_penalty)
    , translated_reference_(NULL)
    , reference_length_(0)
  {
    BuildDefaultMatrix();
  }

  Aligner::Aligner(const int8_t* score_matrix,
                   const int&    score_matrix_size,
                   const int8_t* translation_matrix,
                   const int&    translation_matrix_size)

    : score_matrix_(NULL)
    , score_matrix_size_(score_matrix_size)
    , translation_matrix_(NULL)
    , match_score_(2)
    , mismatch_penalty_(2)
    , gap_opening_penalty_(3)
    , gap_extending_penalty_(1)
    , translated_reference_(NULL)
    , reference_length_(0)
  {
    score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
    memcpy(score_matrix_, score_matrix, sizeof(int8_t) * score_matrix_size_ * score_matrix_size_);
    translation_matrix_ = new int8_t[translation_matrix_size];
    memcpy(translation_matrix_, translation_matrix, sizeof(int8_t) * translation_matrix_size);
  }


  Aligner::Aligner(const Aligner &aligner) {
    match_score_ = aligner.match_score_;
    mismatch_penalty_ = aligner.mismatch_penalty_;
    gap_opening_penalty_ = aligner.gap_opening_penalty_;
    gap_extending_penalty_ = aligner.gap_extending_penalty_;
    ambiguity_penalty_ = aligner.ambiguity_penalty_;

    translated_reference_ = nullptr;
    reference_length_ = 0;

    score_matrix_size_ = aligner.score_matrix_size_;
    score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
    memcpy(score_matrix_, aligner.score_matrix_, sizeof(int8_t) * score_matrix_size_ * score_matrix_size_);

    auto translation_matrix_size = SizeOfArray(kBaseTranslation);
    translation_matrix_ = new int8_t[translation_matrix_size];
    memcpy(translation_matrix_, aligner.translation_matrix_, sizeof(int8_t) * translation_matrix_size);
  }

  Aligner::~Aligner(void){
    Clear();
  }

  int Aligner::SetReferenceSequence(const char* seq, const int& length) {

    int len = 0;
    if (translation_matrix_) {
      // calculate the valid length
      //int calculated_ref_length = static_cast<int>(strlen(seq));
      //int valid_length = (calculated_ref_length > length)
      //                   ? length : calculated_ref_length;
      int valid_length = length;
      // delete the current buffer
      CleanReferenceSequence();
      // allocate a new buffer
      translated_reference_ = new int8_t[valid_length];

      len = TranslateBase(seq, valid_length, translated_reference_);
    } else {
      // nothing
    }

    reference_length_ = len;
    return len;


  }

  int Aligner::TranslateBase(const char* bases, const int& length,
                             int8_t* translated) const {

    const char* ptr = bases;
    int len = 0;
    for (int i = 0; i < length; ++i) {
      translated[i] = translation_matrix_[(int) *ptr];
      ++ptr;
      ++len;
    }

    return len;
  }


  bool Aligner::Align(const char* query, const Filter& filter,
                      Alignment* alignment, const int32_t maskLen) const
  {
    if (!translation_matrix_) return false;
    if (reference_length_ == 0) return false;

    int query_len = strlen(query);
    if (query_len == 0) return false;
    int8_t* translated_query = new int8_t[query_len];
    TranslateBase(query, query_len, translated_query);

    const int8_t score_size = 2;
    s_profile* profile = ssw_init(translated_query, query_len, score_matrix_,
                                  score_matrix_size_, score_size);

    uint8_t flag = 0;
    SetFlag(filter, &flag);
    s_align* s_al = ssw_align(profile, translated_reference_, reference_length_,
                              static_cast<int>(gap_opening_penalty_),
                              static_cast<int>(gap_extending_penalty_),
                              flag, filter.score_filter, filter.distance_filter, maskLen);

    alignment->Clear();
    ConvertAlignment(*s_al, query_len, alignment);
    alignment->mismatches = CalculateNumberMismatch(alignment, translated_reference_, translated_query, query_len);


    // Free memory
    delete [] translated_query;
    align_destroy(s_al);
    init_destroy(profile);

    return true;
  }

  bool Aligner::Align(const char* query, const char* ref, const int& ref_len,
                      const Filter& filter, Alignment* alignment, const int32_t maskLen) const
  {
    return Align(query, strlen(query), ref, ref_len, filter, alignment, maskLen);
  }

  bool Aligner::Align(const char* query, const int& query_len, const char* ref, const int& ref_len,
                      const Filter& filter, Alignment* alignment, const int32_t maskLen) const
  {
    if (!translation_matrix_) return false;

    if (query_len == 0) return false;
    int8_t* translated_query = new int8_t[query_len];
    TranslateBase(query, query_len, translated_query);

    // calculate the valid length
    int valid_ref_len = ref_len;
    int8_t* translated_ref = new int8_t[valid_ref_len];
    TranslateBase(ref, valid_ref_len, translated_ref);


    const int8_t score_size = 2;
    s_profile* profile = ssw_init(translated_query, query_len, score_matrix_,
                                  score_matrix_size_, score_size);

    uint8_t flag = 0;
    SetFlag(filter, &flag);
    s_align* s_al = ssw_align(profile, translated_ref, valid_ref_len,
                              static_cast<int>(gap_opening_penalty_),
                              static_cast<int>(gap_extending_penalty_),
                              flag, filter.score_filter, filter.distance_filter, maskLen);

    alignment->Clear();
    ConvertAlignment(*s_al, query_len, alignment);
    alignment->mismatches = CalculateNumberMismatch(alignment, translated_ref, translated_query, query_len);

    // Free memory
    delete [] translated_query;
    delete [] translated_ref;
    align_destroy(s_al);
    init_destroy(profile);

    return true;
  }

  void Aligner::Clear(void) {
    ClearMatrices();
    CleanReferenceSequence();
  }

  void Aligner::SetAllDefault(void) {
    score_matrix_size_     = 5;
    match_score_           = 2;
    mismatch_penalty_      = 2;
    gap_opening_penalty_   = 3;
    gap_extending_penalty_ = 1;
    ambiguity_penalty_     = 2;
    reference_length_      = 0;
  }

  bool Aligner::ReBuild(void) {
    if (translation_matrix_) return false;

    SetAllDefault();
    BuildDefaultMatrix();

    return true;
  }

  bool Aligner::ReBuild(
    const uint8_t& match_score,
    const uint8_t& mismatch_penalty,
    const uint8_t& gap_opening_penalty,
    const uint8_t& gap_extending_penalty,
    const uint8_t& ambiguity_penalty) {
    if (translation_matrix_) return false;

    SetAllDefault();

    match_score_           = match_score;
    mismatch_penalty_      = mismatch_penalty;
    gap_opening_penalty_   = gap_opening_penalty;
    gap_extending_penalty_ = gap_extending_penalty;
    ambiguity_penalty_     = ambiguity_penalty;

    BuildDefaultMatrix();

    return true;
  }

  bool Aligner::ReBuild(
    const int8_t* score_matrix,
    const int&    score_matrix_size,
    const int8_t* translation_matrix,
    const int&    translation_matrix_size) {

    score_matrix_size_ = score_matrix_size;
    match_score_ = 2;
    mismatch_penalty_ = 2;
    gap_opening_penalty_ = 3;
    gap_extending_penalty_ = 1;
    reference_length_ = 0;

    ClearMatrices();
    score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
    memcpy(score_matrix_, score_matrix, sizeof(int8_t) * score_matrix_size_ * score_matrix_size_);
    translation_matrix_ = new int8_t[translation_matrix_size];
    memcpy(translation_matrix_, translation_matrix, sizeof(int8_t) * translation_matrix_size);

    return true;
  }

  void Aligner::BuildDefaultMatrix(void) {
    ClearMatrices();
    score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
    BuildSwScoreMatrix(match_score_, mismatch_penalty_, ambiguity_penalty_, score_matrix_);
    translation_matrix_ = new int8_t[SizeOfArray(kBaseTranslation)];
    memcpy(translation_matrix_, kBaseTranslation, sizeof(int8_t) * SizeOfArray(kBaseTranslation));
  }

  void Aligner::ClearMatrices(void) {
    delete [] score_matrix_;
    score_matrix_ = NULL;

    delete [] translation_matrix_;
    translation_matrix_ = NULL;
  }
} // namespace StripedSmithWaterman
