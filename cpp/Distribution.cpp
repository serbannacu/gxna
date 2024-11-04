#include "Distribution.h"

#include <cmath>

namespace gxna {

// Helpers for incomplete beta function calculation
// Use the Lentz algorithm to compute continuous fraction

// Avoid division by zero
inline void checkZero(double x) {
    constexpr double Epsilon = 1e-30;
    if (std::fabs(x) < Epsilon)
        x = Epsilon;
}

// Iterative updates of the C and D continuous fraction coefficients
inline void updateD(double& d, double coeff) {
    d = 1 + coeff * d;
    checkZero(d);
    d = 1 / d;
}

inline double updateCD(double& frac, double& c, double& d, double coeff) {
    c = 1 + coeff / c;
    checkZero(c);
    updateD(d, coeff);
    double prod = c * d;
    frac *= prod;
    return prod;
}

// Compute continuous fraction
static double betaLentz(double a, double b, double x) {
    constexpr int MaxIterations = 500;
    constexpr double MaxError = 1e-8;
    double sum = a + b;
    double q = a;
    double c = 1;
    double d = 1;
    updateD(d, -sum / (a + 1) * x);
    double frac = d;
    for (int m = 1; m <= MaxIterations; ++m) {
        double coeff;
        q += 2;
        coeff = m * (--b) / ((q - 1) * q);
        updateCD(frac, c, d, coeff * x);
        coeff = (++a) * (++sum) / (q * (q + 1));
        double val = updateCD(frac, c, d, -coeff * x) - 1;
        if (std::fabs(val) < MaxError)
            break;
    }
    return frac;
}

// Incomplete beta function
// Needs 0 < a, 0 < b, 0 <= x <= 1
double betainc(double a, double b, double x) {
    if (x == 0)
        return 0;
    if (x == 1)
        return 1;
    double logval = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
        + a * std::log(x) + b * std::log(1 - x);
    double val = std::exp(logval);
    if (x < (a + 1) / (a + b + 2))  // continuous fraction converges fast
        return val * betaLentz(a, b, x) / a;
    else
        return 1 - val * betaLentz(b, a, 1 - x) / b;
}

// Inverse error function
// Overflows within 10^-10 or so of 0 and 1
// Max error has order 10^-3
double erfinv(double x) {
    constexpr double A = 7.1422296;
    constexpr double B = 4.546885;  // 2 * A / PI
    double p = std::log(4.0 * x * (1-x));
    double q = 0.5 * p + B;
    return std::sqrt(std::sqrt(q * q - p * A) - q);
}

constexpr double SQRT2 = 1.414214;

double normCDF(double x) {
    constexpr double SQRT2INV = 0.7071067811865475;
    return 1 - 0.5 * std::erfc(x * SQRT2INV);
}

// Computes quantiles for gaussian: f(0) = -Inf, f(1/2) = 0, f(1) = Inf
double normCDFInv(double x) {
    double val = SQRT2 * erfinv(x);
    return x >= 0.5 ? val : -val;
}

double tCDF(double x, double n) {
    return 1 - 0.5 * betainc(n / 2, 0.5, n / (n + x * x));
}

double fCDF(double x, double n1, double n2) {
    return 1 - betainc(n2 / 2, n1 / 2, n2 / (n2 + n1 * x));
}

// Normal equivalent to F statistic, avoids overflow
// Returns y with P(|normal| < y) = P(F < x)
// Needs x >= 0
double zfCDF(double x, double n1, double n2) {
    double y = n2 / (n2 + n1 * x);
    double bInc = betainc(n2 / 2, n1 / 2, y);
    return SQRT2 * erfinv(bInc / 2);
}

// Normal equivalent to T statistic, avoids overflow
// Returns y with P(|normal| < y) = P(T < x)
// Needs x >= 0
double ztCDF(double x, double n) {
    return x >= 0 ? zfCDF(x * x, 1, n) : -zfCDF(x * x, 1, n);
}

}  // namespace gxna
