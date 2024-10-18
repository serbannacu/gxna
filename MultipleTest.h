#pragma once

#include "RandomPermutation.h"
#include "Statistics.h"
#include "VectorUtil.h"

#include <algorithm>
#include <limits>
#include <vector>

namespace gxna {

// Class for testing multiple hypotheses
// It computes raw and adjusted P-values using resampling (permutation tests)
// The method ScoreCalculator::operator() is called for every new permutation and must return a vector of scores
// For two-sided testing (typical), scores should be postive (e.g. the absolute value of the T statistic)
// The first permutation generated MUST be id

template<class ScoreCalculator>
class MultipleTest {
public:
    MultipleTest(ScoreCalculator& calc, int nObjects) : 
        m_calc(calc),
        m_nObjects(nObjects),
        m_rawP(nObjects, 1.0),
        m_adjP(nObjects, 1.0)
    {}

    double getRawP(int i) const { return m_rawP[i]; }
    double getAdjP(int i) const { return m_adjP[i]; }
    double getRank(int i) const { return m_rank[i]; }

    void test(PermutationGenerator& pg, bool scaled) {
        if (scaled)
            maxTscaled(pg);
        else
            maxT(pg);

        // finish computing P-values
        for (int i = 1; i < m_nObjects; i++) {
            if (m_adjP[i] < m_adjP[i-1]) // enforce monotonicity
                m_adjP[i] = m_adjP[i-1];
        }
        m_rawP /= pg.count();
        m_adjP /= pg.count();
    }
    
    void maxT(PermutationGenerator& pg) {
        while (pg.next()) { // generate new permutation
            auto score = m_calc(pg.get());
            if (pg.count() == 1) // first perm housekeeping
                setRank(m_origT = score);
            else
                updateP(score);
        }
    }

    void maxTscaled(PermutationGenerator& pg) {
        std::vector< std::vector<double> > v_score;
        std::vector<FastDataSet> v_data(m_nObjects);

        while (pg.next()) { // compute all scores for all perms
            auto score = m_calc(pg.get());
            v_score.emplace_back(score);
            for (int i = 0; i < m_nObjects; ++i)
                v_data[i].insert(score[i]);
        }

        for (int i = 0; i < m_nObjects; ++i) { // rescale scores
            double mu = v_data[i].mean();
            double sigma = v_data[i].sd();
            for (auto& score : v_score) {
                score[i] -= mu;
                if (sigma > 0)
                    score[i] /= sigma;
            }
        }

        setRank(m_origT = v_score[0]);
        for (int j = 1; j < v_score.size(); ++j)
            updateP(v_score[j]);
    }

private:
    void setRank(const std::vector<double>& score) {
        m_rank.resize(m_nObjects);
        for (int i = 0; i < m_nObjects; i++)
            m_rank[i] = i;
        std::sort(m_rank.begin(), m_rank.end(), [&](int j, int k) { return score[j] > score[k]; });
    }

    void updateP(const std::vector<double>& score) {
        double currentMax = std::numeric_limits<double>::min();
        for (int i = m_nObjects - 1; i >= 0; i--) { // go thru all objects
            int k = m_rank[i];
            double T = score[k];
            if (currentMax < T)
                currentMax = T;
            if (T >= m_origT[k]) // update raw count
                ++m_rawP[i];
            if (currentMax >= m_origT[k]) // update adjusted count
                ++m_adjP[i];
        }
    }
    
    ScoreCalculator& m_calc;
    int m_nObjects;
    std::vector<double> m_origT;
    std::vector<int> m_rank;
    std::vector<double> m_rawP;
    std::vector<double> m_adjP;
};

} // namespace gxna
