// Common distributions (normal, T, F) and special functions
//
// For a random variable X with distribution Dist:
//     DistCDF(x) is Prob(X <= x)
//     DistCDFInv(x) gives the quantile corresponding to x
//     zDistCDF(x) is x converted into an equivalent gaussian

#pragma once

namespace gxna {

// Normal distribution
double normCDF(double x);
double normCDFInv(double x);

// t distribution
double tCDF(double x, double n);
double ztCDF(double x, double n);

// F distribution
double fCDF(double x, double n1, double n2);
double zfCDF(double x, double n1, double n2);

// Special functions used in empirical Bayes calculations
double digamma(double x);
double trigamma(double x);
double trigammainv(double x);

}  // namespace gxna
