#include "RandomPermutation.h"

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
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

// Apply random permutation to array
static void scramble(int *p, size_t len) {
    for (size_t i = 1; i < len; ++i) {
        size_t j = urand(i+1);  // uniform among 0, 1, ..., i
        auto temp = p[i];  // swap p[i] and p[j]
        p[i] = p[j];
        p[j] = temp;
    }
}

void Permutation::randomize() {
    scramble(&m_v[0], m_v.size());
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
    if (m_count >= m_limit)  // done
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

    std::map<int, std::vector<int> > label2positions;
    for (size_t i = 0; i < n; ++i)
        label2positions[label[i]].emplace_back(i);
    for (auto x : label2positions)
        m_labelData.emplace_back(x.second);
}

// For each label, generate a random permutation of 0, 1, ... k - 1, where
// k is the number of elements that have this label. Use this permutation
// to scramble the element positions among themselves. This will generate
// an invariant permutation.

void InvariantPermutation::update() {
    for (auto& data : m_labelData) {
        data.perm.randomize();
        auto& pos = data.positions;
        for (size_t i = 0; i < pos.size(); ++i)  // scramble positions
            m_perm[pos[i]] = pos[data.perm[i]];
    }
}

}  // namespace gxna
