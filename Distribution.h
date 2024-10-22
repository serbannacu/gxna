// Common distributions (normal, T, F) and special functions
//
// For a random variable X with distribution Dist:
//     DistCDF(x) is Prob(X <= x)
//     DistCDFInv(x) gives the quantile corresponding to x
//     zDistCDF(x) is x converted into an equivalent gaussian

#pragma once

namespace gxna {

// Normal distribution
double normCDF(const double x);
double normCDFInv(const double x);

// t distribution
double tCDF(const double x, const double n);
double ztCDF(const double x, const double n);

// F distribution
double fCDF(const double x, const double n1, const double n2);
double zfCDF(const double x, const double n1, const double n2);

// Special functions used in empirical Bayes calculations
double digamma(const double x);
double trigamma(const double x);
double trigammainv(const double y);

} // namespace gxna
