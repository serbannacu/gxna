// Common distributions (normal, T, F) and special functions
//
// For a random variable X with distribution Dist:
//     DistCDF(x) is Prob(X <= x)
//     DistCDFInv(x) gives the quantile corresponding to x
//     zDistCDF(x) is x converted into an equivalent gaussian

#pragma once

namespace gxna {

double normCDF(const double x);
double normCDFInv(const double x);

double tCDF(const double x, const double n);
double ztCDF(const double x, const double n);

double fCDF(const double x, const double n1, const double n2);
double zfCDF(const double x, const double n1, const double n2);

double digamma(const double x);
double trigamma(const double x);
double trigammainv(const double y);

} // namespace gxna
