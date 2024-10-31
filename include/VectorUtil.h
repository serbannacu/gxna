#pragma once

#include <vector>

namespace gxna {

// Basic vector arithmetic: not optimized, but good enough for our purposes, and avoids using third-party libraries
// No size checks

template<class T>
std::vector<T>& operator+=(std::vector<T>& x, const std::vector<T>& y) {
    for (size_t i = 0; i < x.size(); ++i)
        x[i] += y[i];
    return x;
}

template<class T>
std::vector<T>& operator-=(std::vector<T>& x, const std::vector<T>& y) {
    for (size_t i = 0; i < x.size(); ++i)
        x[i] -= y[i];
    return x;
}

template<class T, class U>
std::vector<T>& operator/=(std::vector<T>& x, const U& y) {
    for (size_t i = 0; i < x.size(); ++i)
        x[i] /= y;
    return x;
}

// Vector input/output

template<typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& vec) {
    T x;
    while (is >> x)
        vec.emplace_back(x);
    return is;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    const char sep = ' ';
    for (auto& x : vec)
        os << x << sep;
    return os;
}

} // namespace gxna
