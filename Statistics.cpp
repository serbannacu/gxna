#include "Distribution.h"
#include "Statistics.h"

namespace gxna {

// T statistic between two data sets, different variances assumed
double tstat(FastDataSet& ds0, FastDataSet& ds1) {
    if (ds0.n < 2 || ds1.n < 2)
        return 0;
    else {
        double mu0 = ds0.sumx / ds0.n;
        double mu1 = ds1.sumx / ds1.n;
        double vardiff = (ds1.sumxx / ds1.n - mu1 * mu1) / ds1.n
            + (ds0.sumxx / ds0.n - mu0 * mu0) / ds0.n;
        return vardiff > 0 ? (mu1 - mu0) / sqrt(vardiff) : 0;
    }
}

// T statistic between two data sets, equal variances assumed
double tstateq(FastDataSet& ds0, FastDataSet& ds1) {
    if (ds0.n < 1 || ds1.n < 1)
        return 0;
    else {
        double mu0 = ds0.sumx / ds0.n;
        double mu1 = ds1.sumx / ds1.n;
        int n = ds0.n + ds1.n;
        double mu = (ds0.sumx + ds1.sumx) / n;
        double var = ((ds0.sumxx + ds1.sumxx) / n - mu * mu) / n;
        return var > 0 ? (mu1 - mu0) / sqrt(var) : 0;
    }
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

double fstat(const double *x, const std::vector<int>& pheno, const int nPheno) {
    std::vector<FastDataSet> v_data(nPheno); // stores stats for each phenotype
    for (auto p : pheno)
        v_data[p].insert(*x++);
    double sum = 0, sumsq = 0, sumrss = 0;
    for (auto& data : v_data) {
        sum += data.sumx;
        sumsq += data.sumxx;
        if (data.n > 1)
            sumrss += (data.sumxx - data.sumx * data.sumx / data.n);
    }
    double rssnull = sumsq - sum * sum / pheno.size();
    double rssdiff = rssnull - sumrss;
    if (sumrss > 0 && nPheno > 1)
        return rssdiff / (nPheno - 1) / sumrss * (pheno.size() - nPheno);
    else
        return 0; // not quite correct, but will prevent overflow
}

void EmpiricalBayes::estimate(double logVarMean, double logVarVar, int nGenes, int df) {
    if (nGenes < 2 || df <= 0)
        return;
    double y = logVarVar * (1 + 1.0 / (nGenes - 1)) - trigamma(df / 2.0);
    double val = logVarMean - digamma(df / 2);
    if (y > 0) {
        m_df = 2 * trigammainv(y);
        m_var = exp(val + digamma(m_df / 2) + log(df / m_df));
    }
    else {
        m_df = 100 * df; // just need something large
        m_var = exp(val + log(df / 2));
    }
    m_ratio = m_df / df;
}

double EmpiricalBayes::shrinkageFactor(double geneVar) const {
    double val = geneVar > 0 ? (1 + m_ratio) / (1 + m_ratio * m_var / geneVar) : 1.0;
    return sqrt(val);
}

} // namespace gxna
