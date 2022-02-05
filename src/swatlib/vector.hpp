#ifndef VECTOR_HELPER_H
#define VECTOR_HELPER_H
#include <iostream>
#include <sstream>
#include <utility>

namespace swatlib {

// from: https://stackoverflow.com/a/54383242
template<class TupType, size_t... I>
inline std::ostream& tuple_print(std::ostream& os,
                          const TupType& _tup, std::index_sequence<I...>)
{
    os << "(";
    (..., (os << (I == 0 ? "" : ", ") << std::get<I>(_tup)));
    os << ")";
    return os;
}

template<class... T>
inline std::ostream& operator<< (std::ostream& os, const std::tuple<T...>& _tup)
{
    return tuple_print(os, _tup, std::make_index_sequence<sizeof...(T)>());
}

template<typename T>
inline std::string printVector(std::vector<T> v) {
    std::stringstream ss;
    ss << "Vector[";
    for (auto& vv : v) {
        ss << std::hex << vv << ",";
    }
    ss << "]";
    return ss.str();
}

template<typename T>
inline std::string printVectorD(std::vector<T> v) {
    std::stringstream ss;
    ss << "Vector[";
    for (auto& vv : v) {
        ss << vv << ",";
    }
    ss << "]";
    return ss.str();
}

}

#endif // VECTOR_HELPER_H