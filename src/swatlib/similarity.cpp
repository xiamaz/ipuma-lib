#include "similarity.h"


int swatlib::simpleSimilarity(char a, char b) {
    if (a == b)
        return 1;
    return -1;
}

void swatlib::to_json(nlohmann::json& j, const swatlib::Similarity& s) {
	j = swatlib::similarityToStr(s);
}

void swatlib::from_json(const nlohmann::json& j, swatlib::Similarity& s) {
	s = swatlib::strToSimilarity(j);
}

template <typename T>
swatlib::Matrix<T> swatlib::scaleMatrix(swatlib::Matrix<T> mat, T pos, T neg, T ambig) {
    auto [m, n] = mat.shape();
    int size = m * n;
    std::vector<T> scaled(size);
    int8_t* mv = mat.data();
    for (int i = 0; i < size; ++i) {
        scaled[i] = mv[i] ? pos : neg;
    }
    // insert ambiguity for N values
    for (int i = 0; i < m; ++i) {
        scaled[i] = ambig;
        scaled[i*n] = ambig;
    }
    return swatlib::Matrix<T>(m, n, scaled);
}

swatlib::Matrix<int8_t> swatlib::createQueryProfile(const std::vector<uint8_t>& sequence, const swatlib::Matrix<int8_t>& scores) {
    auto [m, n] = scores.shape();
    swatlib::Matrix<int8_t> profile(sequence.size(), m);
    for (int i = 0; i < sequence.size(); ++i) {
        for (int j = 0; j < n; ++j) {
            profile(i, j) = scores(sequence[i], j);
        }
    }
    return profile;
}

swatlib::Similarity swatlib::strToSimilarity(std::string s) {
    if (s == "na") {
        return swatlib::Similarity::nucleicAcid;
    } else if (s == "blosum50") {
        return swatlib::Similarity::blosum50;
    } else if (s == "blosum62") {
        return swatlib::Similarity::blosum62;
    } else if (s == "simple") {
        return swatlib::Similarity::simple;
    } else {
        throw std::runtime_error("Unknown similarity matrix type.");
    }
}
std::string swatlib::similarityToStr(swatlib::Similarity s) {
    switch(s) {
        case swatlib::Similarity::nucleicAcid:
        return "na";
        break;
        case swatlib::Similarity::blosum62:
        return "blosum62";
        break;
        case swatlib::Similarity::blosum50:
        return "blosum50";
        break;
        case swatlib::Similarity::simple:
        return "simple";
        break;
    }
    throw std::runtime_error("Invalid similarity type.");
}

template <typename T>
swatlib::Matrix<T> convertMatrix(swatlib::Matrix<int8_t> om) {
    auto [m, n] = m.shape();
    swatlib::Matrix<T> cm(m, n);
    T* cmp = cm.data();
    int8_t* omp = om.data();
    for (int i = 0; i < m * n; ++i) {
        cmp[i] = static_cast<T>(omp[i]);
    }
    return cm;
}

swatlib::Matrix<int8_t> swatlib::selectMatrix(swatlib::Similarity s, int8_t matchValue, int8_t mismatchValue, int8_t ambiguityValue){
    switch(s) {
    case swatlib::Similarity::nucleicAcid:
        return swatlib::scaleMatrix<int8_t>(swatlib::NA_MATRIX.copy(), matchValue, mismatchValue, ambiguityValue);
        break;
    case swatlib::Similarity::blosum50:
        return swatlib::BLOSUM50_MATRIX.copy();
        break;
    case swatlib::Similarity::blosum62:
        return swatlib::BLOSUM62_MATRIX.copy();
        break;
    case swatlib::Similarity::simple:
        throw std::runtime_error("No matrix for simple penalty available.");
        break;
    }
    throw std::runtime_error("Invalid similarity type.");
}