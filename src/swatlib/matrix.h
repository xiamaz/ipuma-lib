#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>
#include <iostream>
#include <tuple>

namespace swatlib {

template<typename T>
class Matrix {
protected:
    std::vector<T> vdata;
    size_t m, n;
public:
    Matrix();
    Matrix(size_t m, size_t n);
    Matrix(size_t m, size_t n, T v);
    Matrix(size_t m, size_t n, std::vector<T> v);

    std::tuple<size_t, size_t> shape() const;
    size_t getM() const;
    size_t getN() const;

    T& get(int i, int j);

    void set(int i, int j, T v);

    /**
     * Return the index of the largest element in the matrix.
     * 
     * If multiple elements are the largest, any of them might be returned.
     */
    std::tuple<int, int> argmax() const;

    /**
     * Copy values from the given vector into own data.
     * 
     * Throws if the lengths of the vectors do not match.
     */
    void copyValues(std::vector<T> newData);

    Matrix copy() const;

    T* data();
    const T* data() const;

    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();

    template<class Y>
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

    /**
     * Get/Set without bounds check.
     */
    T& operator()(int i, int j);
    T operator()(int i, int j) const;

    const std::vector<T> vector() const;

    std::string toString();
    std::vector<std::vector<T>> toVector();
};

}

#endif // MATRIX_H