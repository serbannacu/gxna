#include "RandomPermutation.h"

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>

namespace gxna {

static std::mt19937 mt;  // Mersenne Twister random number generator

inline int urand(int n) {  // uniform random between 0 and n-1
    static std::uniform_real_distribution<> dist(0.0, 1.0);
    return dist(mt) * n;
}

void Permutation::seed(int val) {
    mt.seed(val);
}

void Permutation::randomize(size_t start, size_t len) {
    auto p = &m_v[start];
    for (size_t i = 1; i < len; ++i) {
        size_t j = urand(i+1);  // uniform among 0, 1, ..., i
        auto temp = p[i];  // swap p[i] and p[j]
        p[i] = p[j];
        p[j] = temp;
    }
}

void PermutationHistogram::insert(const Permutation& p) {
    assert(p.size() == m_n);
    for (size_t i = 0; i < m_n; ++i)
        ++m_count[i][p[i]];
    ++m_n_perm;
}

void PermutationHistogram::print(std::ostream& os, int prec) const {
    os << std::fixed << std::setprecision(prec);
    for (auto& vec : m_count) {
        for (auto& val : vec)
            os << val / double(m_n_perm) << ' ';
        os << '\n';
    }
}

bool PermutationGenerator::next() {
    if (m_count >= m_limit)
        return false;
    if (m_count++)  // first next() call leaves Id permutation in place
        update();
    if (m_verbose)
        showProgress();
    return true;
}

void PermutationGenerator::showProgress() {
    int val = 100 * m_count / m_limit;
    if (m_percentage != val) {
        m_percentage = val;
        std::cerr << "Permutation " << m_count << '/' << m_limit << ' '
                  << val << '%' << (m_count == m_limit ? '\n' : '\r');
    }
}

InvariantPermutation::InvariantPermutation(size_t n, size_t limit,
                                           const std::vector<int>& label)
    : PermutationGenerator(n, limit) {
    assert(n == label.size());
    int k = -1;  // current label
    for (auto val : label) {
        if (m_labelCount.empty() || k != val) {
            k = val;
            m_labelCount.emplace_back(0);
        }
        ++m_labelCount.back();
    }
}

void InvariantPermutation::update() {
    size_t sum = 0;
    for (auto val : m_labelCount) {
        m_perm.randomize(sum, val);
        sum += val;
    }
}

}  // namespace gxna
