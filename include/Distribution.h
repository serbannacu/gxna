#pragma once

namespace gxna {

// Common distributions (normal, T, F)

// For a random variable X with distribution dist:
//     distCDF(x) is Prob(X <= x)
//     distCDFInv(x) gives the quantile corresponding to x
//     zdistCDF(x) is x converted into an equivalent gaussian

double betainc(double a, double b, double x); // incomplete beta function
double erfinv(double x); // inverse error function

// Normal distribution
double normCDF(double x);
double normCDFInv(double x);

// t distribution
double tCDF(double x, double n);
double ztCDF(double x, double n);

// F distribution
double fCDF(double x, double n1, double n2);
double zfCDF(double x, double n1, double n2);

}  // namespace gxna
