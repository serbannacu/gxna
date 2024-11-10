#pragma once

namespace gxna {

// Common distributions: normal, t, F.

// Regularized incomplete beta function (beta CDF).
// Needs 0 < a, 0 < b, 0 <= x <= 1.

double betainc(double a, double b, double x);

// Cumulative Distribution and Survival Functions.
// CDF(x) is Prob(X <= x).
// SF(x) is Prob(X > x) == 1 - CDF(x).
// For tail probabilities, CDF is more accurate for small x abd SF for large x.

double zCDF(double z);  // Gaussian
double zSF(double z);

double tCDF(double t, double n);
double tSF(double t, double n);

double fCDF(double f, double n1, double n2);
double fSF(double f, double n1, double n2);

// Inverse CDF; need 0 <= y <= 1.

double zCDFinv(double y);

// Convert t- or F- score into z-score.
// In the formulas below z is the function return value and Z is a standard normal.

double t2z(double t, double n);  // P(T < t) == P(Z < z)
double f2z(double f, double n1, double n2);  // P(F < f) == P(|Z| < z); needs f >= 0

}  // namespace gxna
