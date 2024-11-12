#include "RandomPermutation.h"
#include "VectorUtil.h"

#include <cassert>
#include <cmath>
#include <utility>  // make_pair
#include <vector>
#include "gtest/gtest.h"

namespace gxna {

// Given a vector of observed frequencies from sampling a discrete uniform
// distribution, check whether they are close enough to the expected value.

template<typename T>
T checkDiff(const std::vector<T>& v, double expected, double limit) {
    T maxDiff = 0;
    for (T observed : v) {
        auto val = std::abs(observed - expected);
        if (maxDiff < val)
            maxDiff = val;
    }
    EXPECT_LE(maxDiff, limit) << ::testing::PrintToString(v) << ' ' << expected;
    return maxDiff;
}

// Mean and SE of observed frequency from sampling a discrete uniform distribution.

auto getMeanAndSE(double nValues, double nSamples) {
    auto mean = nSamples / nValues;  // expected count
    auto var = mean * (1 - 1 / nValues);
    auto se = std::sqrt(var);  // standard error
    return std::make_pair(mean, se);
}

void testUniform(size_t nValues, size_t nSamples, double maxRatio) {
    PermutationHistogram hist(nValues);
    UniformPermutation pg(nValues, nSamples);
    while (pg.next())
        hist.insert(pg.get());

    // For position i of random perm p, check p[i] is uniform.

    auto [mean, se] = getMeanAndSE(nValues, nSamples);
    for (size_t i = 0; i < nValues; ++i)
        checkDiff(hist.getCount(i), mean, maxRatio * se);
}

void testInvariant(size_t nValues, size_t nSamples, double maxRatio, std::vector<int> label) {
    if (label.empty())
        label.resize(nValues);
    else
        assert(label.size() == nValues);

    PermutationHistogram hist(nValues);
    InvariantPermutation pg(nValues, nSamples, label);
    while (pg.next())
        hist.insert(pg.get());

    // For position i of random perm p, check p[i] is uniform over the positions
    // that have the same label.

    for (size_t i = 0; i < nValues; ++i) {
        auto count = hist.getCount(i);
        std::vector<double> v;
        for (size_t j = 0; j < nValues; ++j) {
            if (label[j] == label[i])
                v.emplace_back(count[j]);
            else  // invariant perm must preserve label
                EXPECT_EQ(count[j], 0);
        }
        auto [mean, se] = getMeanAndSE(v.size(), nSamples);
        checkDiff(v, mean, maxRatio * se);
    }
}

TEST(Permutation, Uniform) {
    testUniform(10, 1000000, 3.0);
    testUniform(15, 100000, 2.8);
}

TEST(Permutation, Invariant) {
    testInvariant(10, 1000000, 3.0, {});
    testInvariant(10, 100000, 3.0, {5, 5, 5, 1, 1, 1, 1, 2, 2, 2});
    testInvariant(10, 100000, 3.0, {2, 5, 5, 1, 1, 1, 1, 0, 2, 2});
}

}  // namespace gxna
