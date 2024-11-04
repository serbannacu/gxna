#pragma once

#include <cmath>
#include <vector>

namespace gxna {

// Statistical tools: mean, variance, T-stat, F-stat.

class FastDataSet {
 public:
    FastDataSet() {}

    FastDataSet(const std::vector<double>& vec) {
        insert(vec);
    }

    void insert(double x) {  // add data point
        m_x += x;
        m_xx += x * x;
        ++m_n;
    }

    void insert(const std::vector<double>& vec) {
        for (double x : vec)
            insert(x);
    }

    FastDataSet& operator+=(const FastDataSet& rhs) {
        m_x += rhs.m_x;
        m_xx += rhs.m_xx;
        m_n += rhs.m_n;
        return *this;
    }

    int size() const {
        return m_n;
    }

    // Safe methods (with size checks)

    double mean() const {
        return m_n ? m_x / m_n : 0;
    }

    double var() const {
        return m_n > 1 ? getVar() : 0;
    }

    double sd() const {
        return m_n > 1 ? std::sqrt(getVar()) : 0;
    }

    double se() const {  // standard error
        return m_n > 1 ? std::sqrt(getVar() / m_n) : 0;
    }

    double rss() const {  // residual sum of squares; useful for t and F tests
        return m_n > 1 ? m_xx - m_x * m_x / m_n : 0;
    }

    friend FastDataSet operator+(FastDataSet x, const FastDataSet& y) {
        return x += y;
    }

    // Unsafe methods (no size checks so a bit faster)

    void meanVar(double& mean_, double& var_) const {
        mean_ = m_x / m_n;
        var_ = getVar(mean_);
    }

 private:
    double getVar(double mean_) const {  // note the var estimate is the MLE (biased)
        double retval = m_xx / m_n - mean_ * mean_;
        return retval > 0 ? retval : 0;
    }

    double getVar() const {
        return getVar(m_x / m_n);
    }

    int m_n = 0;
    double m_x = 0; // sum of elements
    double m_xx = 0;  // sum of squared elements
};

// t and F statistics
double tstat(const FastDataSet&, const FastDataSet&);  // unequal variances (Welch)
double tstatEqual(const FastDataSet&, const FastDataSet&);  // equal variances (Student)
double fstat(const std::vector<FastDataSet>&);

// t and F statistics for an array with associated phenotypes

// phenotypes 0 and 1 define the two data sets, everything else is ignored
double tstatPheno(const double *x, const std::vector<int>& pheno, int pheno0, int pheno1);

// phenotypes are required to be integers btwn 0 and nPheno - 1
double fstatPheno(const double *x, const std::vector<int>& pheno, const int nPheno);


// Special functions used in empirical Bayes calculations
// All expect x > 0

double harmonic(int n);  // approximate sum(1 / k) from 1 to n
double harmonic2(int n);  // approximate sum(1 / k^2) from 1 to n
double digamma(double x);
double trigamma(double x);

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
