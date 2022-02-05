#include "matrix.h"
#include "vector.hpp"
#include <cstdint>
#include <sstream>
#include <algorithm>

template<typename T>
swatlib::Matrix<T>::Matrix() : vdata(0), m{0}, n{0} {}

template<typename T>
swatlib::Matrix<T>::Matrix(size_t m, size_t n) : vdata(m * n), m{m}, n{n} {}

template<typename T>
swatlib::Matrix<T>::Matrix(size_t m, size_t n, T v) : vdata(m * n, v), m{m}, n{n} {}

template<typename T>
swatlib::Matrix<T>::Matrix(size_t m, size_t n, const std::vector<T> v) : vdata(v), m{m}, n{n} {
    if (m * n != v.size()) {
        throw std::runtime_error("Dimensions do not match. m: " + std::to_string(m) + " n: " + std::to_string(n) + " v: " + std::to_string(v.size()));
    }
}

template<typename T>
std::tuple<size_t, size_t> swatlib::Matrix<T>::shape() const {return {m, n};}

template<typename T>
size_t swatlib::Matrix<T>::getM() const { return m; }

template<typename T>
size_t swatlib::Matrix<T>::getN() const { return n; }


template<typename T>
T& swatlib::Matrix<T>::get(int i, int j) {
    if (i >= m || i < 0 || j >= n || j < 0)
        throw std::out_of_range("Matrix get error.");
    // return data.at(i * n + j);
    return vdata[i * n + j];
}

template<typename T>
void swatlib::Matrix<T>::set(int i, int j, T v) {
    if (i >= m || i < 0 || j >= n || j < 0)
        throw std::out_of_range("Matrix get error.");
    // return data.at(i * n + j);
    vdata[i * n + j] = v;
}

template<typename T>
std::tuple<int, int> swatlib::Matrix<T>::argmax() const {
    int dataIndex = std::distance(vdata.begin(), std::max_element(vdata.begin(), vdata.end()));

    // convert data index into our current shape
    int row = dataIndex / n;
    int col = dataIndex % n;
    return {row, col};
}

template<typename T>
void swatlib::Matrix<T>::copyValues(std::vector<T> newData) {
    if (newData.size() != vdata.size())
        throw std::invalid_argument("Size of data does not match.");

    vdata = newData;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const swatlib::Matrix<T>& m) {
    return os << "Matrix(" << m.m << "," << m.n << ")";
}

template<typename T>
T& swatlib::Matrix<T>::operator()(int i, int j) {
    return vdata[i * n + j];
}

template<typename T>
T swatlib::Matrix<T>::operator()(int i, int j) const {
    return vdata[i * n + j];
}

template<typename T>
typename std::vector<T>::iterator swatlib::Matrix<T>::begin() {
    return vdata.begin();
}

template<typename T>
typename std::vector<T>::iterator swatlib::Matrix<T>::end() {
    return vdata.end();
}

template<typename T>
T* swatlib::Matrix<T>::data() {
    return vdata.data();
}

template<typename T>
const std::vector<T> swatlib::Matrix<T>::vector() const {
    return vdata;
}

template<typename T>
const T* swatlib::Matrix<T>::data() const {
    return vdata.data();
}

template<typename T>
swatlib::Matrix<T> swatlib::Matrix<T>::copy() const {
    return swatlib::Matrix(m, n, vdata);
}

template<typename T>
std::string swatlib::Matrix<T>::toString() {
    std::stringstream ss;
    ss << "Matrix[\n";
    for (int i = 0; i < m; i++) {
        ss << "[";
        for (int j = 0; j < n; j++) {
            if (j > 0)
                ss << ",";
            ss << get(i, j);
        }
        ss << "]\n";
    }
    ss << "]";
    return ss.str();
}

template<typename T>
std::vector<std::vector<T>> swatlib::Matrix<T>::toVector() {
    std::vector<std::vector<T>> v(m);
    for (int i = 0; i < m; ++i) {
        v[i] = std::vector<T>(n);
        for (int j = 0; j < n; ++j) {
            v[i][j] = operator()(m, n);
        }
    }
    return v;
}

template class swatlib::Matrix<int>;
template class swatlib::Matrix<short>;
template class swatlib::Matrix<char>;
template class swatlib::Matrix<int8_t>;
template class swatlib::Matrix<uint8_t>;