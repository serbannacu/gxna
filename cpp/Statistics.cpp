#include "Distribution.h"
#include "Statistics.h"

namespace gxna {

static double safe_div_sqrt(double x, double y) {
    if (y > 0)
        return x / std::sqrt(y);
    else
        return 0;
}

double tstat(const FastDataSet& d1, const FastDataSet& d2) {
    int n1 = d1.size();
    int n2 = d2.size();
    if (n1 < 2 || n2 < 2)
        return 0;
    double mean1, mean2, var1, var2;
    d1.meanVar(mean1, var1);
    d2.meanVar(mean2, var2);
    double val = var1 / (n1 - 1) + var2 / (n2 - 1);
    return safe_div_sqrt(mean2 - mean1, val);
}

double tstatEqual(const FastDataSet& d1, const FastDataSet& d2) {
    int n1 = d1.size();
    int n2 = d2.size();
    int n = n1 + n2;
    if (n1 == 0 || n2 == 0 || n == 2)
        return 0;
    double val = (d1.rss() + d2.rss()) / (n - 2) * (1.0 / n1 + 1.0 / n2);
    return safe_div_sqrt(d2.mean() - d1.mean(), val);
}

double fstat(const std::vector<FastDataSet>& v_data) {
    double rss = 0;
    FastDataSet all;
    for (auto& data : v_data) {
        all += data;
        rss += data.rss();
    }
    auto n = v_data.size();
    if (rss > 0 && n > 1)
        return (all.rss() / rss - 1) * (all.size() - n) / (n - 1);
    else
        return 0;  // not quite correct, but will prevent overflow
}

double tstatPheno(const double *x, const std::vector<int>& pheno, int pheno0, int pheno1) {
    FastDataSet ds0, ds1;
    for (auto p : pheno) {
        if (p == pheno0)
            ds0.insert(*x);
        else if (p == pheno1)
            ds1.insert(*x);
        ++x;
    }
    return tstat(ds0, ds1);
}

double fstatPheno(const double *x, const std::vector<int>& pheno, const int nPheno) {
    std::vector<FastDataSet> v_data(nPheno);
    for (auto p : pheno)
        v_data[p].insert(*x++);
    return fstat(v_data);
}

[[maybe_unused]]
double harmonic(int n) {
    constexpr double Gamma = 0.577215664901532;  // Euler's constant
    double inv = 1.0 / n;
    double inv2 = inv * inv;
    double inv4 = inv2 * inv2;
    return std::log(n) + Gamma + 0.5 * inv - inv2 / 12 + inv4 / 120;
}

[[maybe_unused]]
double harmonic2(int n) {
    constexpr double Zeta2 = 1.6449340668482264;  // PI ^ 2 / 6
    double inv = 1.0 / n;
    double inv2 = inv * inv;
    double inv4 = inv2 * inv2;
    return Zeta2 - inv * (1 - 0.5 * inv + inv2 / 6 - inv4 / 30);
}

// Gamma special functions are not optimized, but the implementations are fast enough for our use case

double digamma(double x) {
    double sum = 0;
    while (x < 20)
        sum += 1 / x++;
    double inv = 1.0 / x;
    double inv2 = inv * inv;
    double inv4 = inv2 * inv2;
    double val = std::log(x) - 0.5 * inv - inv2 / 12 + inv4 / 120;
    return val - sum;
}

double trigamma(double x) {
    double sum = 0;
    while (x < 20) {
        sum += 1 / (x * x);
        ++x;
    }
    double inv = 1.0 / x;
    double inv2 = inv * inv;
    double inv4 = inv2 * inv2;
    double val = inv * (1 + 0.5 * inv + inv2 / 6 - inv4 / 30);
    return val + sum;
}

// Trigammainv implementation uses bisection. Newton iteration would be faster.
double trigammainv(double y) {
    // Use asymptotics to avoid overflow or underflow
    if (y < 1e-6)
        return 1 / y + 0.5;
    if (y > 1e6)
        return 1 / std::sqrt(y);
    double x = 0.5 + 1 / y;
    double lo = x, hi = x;
    while (trigamma(lo) < y)
        lo /= 2;
    while (trigamma(hi) > y)
        hi *= 2;
    constexpr double eps = 1e-6;
    while (hi - lo > eps) {
        x = 0.5 * (lo + hi);
        double val = trigamma(x);
        if (val == y)
            break;
        else if (val > y)
            lo = x;
        else
            hi = x;
    }
    return 0.5 * (lo + hi);
}

void EmpiricalBayes::estimate(double logVarMean, double logVarVar, int nGenes, int df) {
    if (nGenes < 2 || df <= 0)
        return;
    double y = logVarVar * (1 + 1.0 / (nGenes - 1)) - trigamma(df / 2.0);
    double val = logVarMean - digamma(df / 2);
    if (y > 0) {
        m_df = 2 * trigammainv(y);
        m_var = std::exp(val + digamma(m_df / 2) + std::log(df / m_df));
    }
    else {
        m_df = 100 * df;  // just need something large
        m_var = std::exp(val + std::log(df / 2));
    }
    m_ratio = m_df / df;
}

double EmpiricalBayes::shrinkageFactor(double geneVar) const {
    if (geneVar > 0) {
        double val = (1 + m_ratio) / (1 + m_ratio * m_var / geneVar);
        return std::sqrt(val);
    }
    else
        return 1.0;
}

}  // namespace gxna
