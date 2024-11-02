#pragma once

#include <cmath>
#include <vector>

namespace gxna {

// Statistical tools: mean, variance, T-stat, F-stat.

class FastDataSet {
 public:
    FastDataSet() {}

    FastDataSet(const std::vector<double>& vec) {
        for (double x : vec)
            insert(x);
    }

    void insert(double x) {  // add data point
        sumx += x;
        sumxx += x * x;
        ++n;
    }

    double mean() const {
        return n > 0 ? sumx / n : 0;
    }

    double var() const {  // note the var estimate is the MLE (biased)
        if (n <= 1) return 0;
        double val = sumx / n;
        return sumxx / n - val * val;
    }

    double sd() const {
        if (n <= 1) return 0;
        double val = var();
        return val > 0 ? sqrt(val) : 0;
    }

    double se() const {
        if (n <= 1) return 0;
        double val = var();
        return val > 0 ? sqrt(val / n) : 0;
    }

    int size() {
        return n;
    }

    friend double tstat(const FastDataSet&, const FastDataSet&);
    friend double tstatEqual(const FastDataSet&, const FastDataSet&);
    friend double fstat(const double *x, const std::vector<int>& pheno, const int nPheno);

 private:
    int n = 0;
    double sumx = 0;
    double sumxx = 0;  // sum of squares
};

// T statistic between two data sets
double tstat(const FastDataSet&, const FastDataSet&);  // unequal variances
double tstatEqual(const FastDataSet&, const FastDataSet&);  // equal variances

// T statistic for a vector with associated phenotypes
// phenotypes 0 and 1 define the two data sets, everything else is ignored
double tstatPheno(const double *x, const std::vector<int>& pheno, int pheno0, int pheno1);

// F statistic for an array with associated phenotypes
// phenotypes are required to be integers btwn 0 and nPheno - 1
double fstat(const double *x, const std::vector<int>& pheno, const int nPheno);

// Empirical Bayes shrinkage based on the following paper:
// Gordon Smyth (2004), Linear models and empirical Bayes methods
// for assessing differential expression in microarray experiments
class EmpiricalBayes {
 public:
    // df is data degrees of freedom
    void estimate(double logVarMean, double logVarVar, int nGenes, int df);

    double shrinkageFactor(double geneVar) const;
    double df() const { return m_df; }
    double var() const { return m_var; }

 private:
    double m_df = 0;
    double m_var = 0;
    double m_ratio = 0;
};

}  // namespace gxna
