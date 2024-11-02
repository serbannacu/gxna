#include "RandomPermutation.h"

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace gxna {

inline int urand(int n) {  // uniform random between 0 and n-1
    return (std::rand() / (((double) RAND_MAX) + 1)) * n;
}

void Permutation::randomize(size_t start, size_t len) {
    auto p = &m_v[start];
    for (size_t i = 1; i < len; ++i) {
        size_t j = urand(i+1);  // uniform among 0, 1, ..., i
        auto temp = p[i];
        p[i] = p[j];
        p[j] = temp;
    }
}

void PermutationHistogram::print(std::ostream& os, int prec) const {
    os << std::fixed << std::setprecision(prec);
    for (size_t i = 0; i < m_n; ++i) {
        for (size_t j = 0; j < m_n; ++j)
            os << m_count[i][j] / double(m_n_add) << ' ';
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

InvariantPermutation::InvariantPermutation(size_t n, size_t limit, const std::vector<int>& types)
    : PermutationGenerator(n, limit) {
    assert(n == types.size());
    int type = 0;
    for (auto val : types) {
        if (m_typeCount.empty() || type != val) {
            type = val;
            m_typeCount.emplace_back(0);
        }
        ++m_typeCount.back();
    }
}

void InvariantPermutation::update() {
    size_t sum = 0;
    for (auto val : m_typeCount) {
        m_perm.randomize(sum, val);
        sum += val;
    }
}

}  // namespace gxna
