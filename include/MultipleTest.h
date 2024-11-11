#pragma once

#include "RandomPermutation.h"
#include "Statistics.h"
#include "VectorUtil.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

namespace gxna {

// Class for testing multiple hypotheses.
// It computes raw and adjusted P-values using resampling (permutation tests).

// The method ScoreCalculator::operator() is called for every new permutation
// and must return a vector of scores.
// Scores should be non-negative, so caller should take the absolute value
// when using two-sided statistics such as T.
// The first permutation generated MUST be id.

// In some use cases, floating point errors in score calculation can result in
// inconsistent behavior across platforms. To prevent this, the ctor argument
// eps can be set to a value that is small but exceeds the expected error.

template<class ScoreCalculator>
class MultipleTest {
 public:
    MultipleTest(ScoreCalculator& calc, int nObjects, double eps = 0)
        : m_calc(calc),
          m_nObjects(nObjects),
          m_eps(eps),
          m_rawP(nObjects, 1.0),
          m_adjP(nObjects, 1.0)
    {}

    double getRawP(int i) const { return m_rawP[i]; }
    double getAdjP(int i) const { return m_adjP[i]; }
    double getRank(int i) const { return m_rank[i]; }

    // Main entry point
    void test(PermutationGenerator& pg, bool scaled) {
        // Compute raw and adjusted counts.
        if (scaled)
            maxTscaled(pg);
        else
            maxT(pg);

        // Enforce monotonicity for adjusted counts.
        for (size_t i = 1; i < m_adjP.size(); ++i) {
            if (m_adjP[i] < m_adjP[i-1])
                m_adjP[i] = m_adjP[i-1];
        }

        // Convert counts into P-values.
        m_rawP /= pg.count();
        m_adjP /= pg.count();
    }

    void maxT(PermutationGenerator& pg) {
        while (pg.next()) {  // generate new permutation
            auto score = m_calc(pg.get());
            assert(score.size() == m_nObjects);
            if (pg.count() == 1)  // first perm housekeeping
                initRank(score);
            else
                updateP(score);
        }
    }

    void maxTscaled(PermutationGenerator& pg) {
        std::vector< std::vector<double> > v_score;
        std::vector<FastDataSet> v_data(m_nObjects);  // to compute mean and sd

        // Compute all scores for all permutations.
        while (pg.next()) {
            auto score = m_calc(pg.get());
            assert(score.size() == m_nObjects);
            v_score.emplace_back(score);
            for (size_t i = 0; i < v_data.size(); ++i)
                v_data[i].insert(score[i]);
        }

        // Rescale each object score by its mean and sd across samples.
        for (size_t i = 0; i < v_data.size(); ++i) {
            double mu = v_data[i].mean();
            double sigma = v_data[i].sd();
            for (auto& score : v_score) {
                score[i] -= mu;
                if (sigma > 0)
                    score[i] /= sigma;
            }
        }

        initRank(v_score[0]);
        for (size_t j = 1; j < v_score.size(); ++j)
            updateP(v_score[j]);
    }

 private:

    // Rank objects in decreasing order of their score.
    // Use stable sort to get consistent behavior across platforms.

    void initRank(const std::vector<double>& score) {
        m_origT = score;
        m_rank.resize(score.size());
        for (size_t i = 0; i < m_rank.size(); ++i)
            m_rank[i] = i;
        std::stable_sort(m_rank.begin(), m_rank.end(),
                         [&](int j, int k) { return score[j] > score[k] + m_eps; });
    }

    // Iterate over objects in increasing order of their original score.
    // Compare resampled scores with original scores, update raw and adjusted
    // counts as needed.

    void updateP(const std::vector<double>& score) {
        double currentMax = std::numeric_limits<double>::min();
        for (int i = m_rank.size() - 1; i >= 0; --i) {  // go thru all objects
            int k = m_rank[i];
            double T = score[k];
            if (currentMax < T)
                currentMax = T;
            if (T >= m_origT[k])  // update raw count
                ++m_rawP[i];
            if (currentMax >= m_origT[k])  // update adjusted count
                ++m_adjP[i];
        }
    }

    ScoreCalculator& m_calc;
    size_t m_nObjects;
    double m_eps;
    std::vector<double> m_origT;
    std::vector<int> m_rank;
    std::vector<double> m_rawP;
    std::vector<double> m_adjP;
};

}  // namespace gxna
