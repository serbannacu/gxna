#include "RandomPermutation.h"
#include "VectorUtil.h"

#include "gtest/gtest.h"

#include <cassert>
#include <cmath>
#include <map>
#include <utility>  // make_pair

namespace gxna {

template<typename T>
T maxabs(const std::vector<T>& v) {
    T best = 0;
    for (T val : v) {
        val = std::abs(val);
        if (best < val)
            best = val;
    }
    return best;
}

auto getMeanAndSE(size_t k, size_t n) {
    auto mean = double(n) / k;  // expected count
    auto var = double(n) * (k - 1) / (k * k);
    auto se = std::sqrt(var);  // standard error
    return std::make_pair(mean, se);
}

// k is permutation size, n is number of permutations

void testUniform(size_t k, size_t n, double maxRatio) {
    PermutationHistogram hist(k);
    UniformPermutation pg(k, n);
    while (pg.next())
        hist.insert(pg.get());

    auto [mean, se] = getMeanAndSE(k, n);
    for (size_t i = 0; i < k; ++i) {
        std::vector<double> v(k, -mean);  // expected
        v += hist.getCount(i);  // observed - expected
        EXPECT_LT(maxabs(v), maxRatio * se);
    }
}

void testInvariant(size_t k, size_t n, double maxRatio, std::vector<int> label) {
    if (label.empty())
        label.resize(k);
    else
        assert(label.size() == k);
    std::map<int, int> label2count;
    for (int val : label)
        ++label2count[val];

    PermutationHistogram hist(k);
    InvariantPermutation pg(k, n, label);
    while (pg.next())
        hist.insert(pg.get());

    for (size_t i = 0; i < k; ++i) {
        auto count = hist.getCount(i);
        std::vector<double> v;
        for (size_t j = 0; j < k; ++j) {
            if (count[j]) {  // must preserve label
                EXPECT_EQ(label[i], label[j]);
                v.emplace_back(count[j]);
            }
        }
        int m = label2count[label[i]];
        EXPECT_EQ(v.size(), m);
    }
}

TEST(Permutation, Uniform) {
    testUniform(10, 1000000, 3.0);
    testUniform(15, 100000, 4.0);
}

TEST(Permutation, Invariant) {
    testInvariant(10, 100000, 3.0, {});
    testInvariant(10, 100000, 3.0, {5, 5, 5, 1, 1, 1, 1, 2, 2, 2});
    testInvariant(10, 100000, 3.0, {2, 5, 5, 1, 1, 1, 1, 2, 2, 2});
}

}  // namespace gxna
