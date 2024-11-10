#include "Distribution.h"

#include <cmath>

namespace gxna {

// Helpers for incomplete beta function calculation.

// Avoid division by zero.

inline void checkZero(double x) {
    constexpr double Epsilon = 1e-30;
    if (std::fabs(x) < Epsilon)
        x = Epsilon;
}

// Compute continuous fraction using the Lentz algorithm.
// The coefficients C_i and D_i are updated iteratively.
// After n iterations, the fraction is estimated as
//   frac = (C_1 D_1) * (C_2 D_2) * ... * (C_n D_n)

inline void updateD(double& d, double coeff) {
    d = 1 + coeff * d;
    checkZero(d);
    d = 1 / d;
}

inline double updateFraction(double& frac, double& c, double& d, double coeff) {
    c = 1 + coeff / c;
    checkZero(c);
    updateD(d, coeff);
    double prod = c * d;
    frac *= prod;
    return prod;
}

// Compute continuous fraction.

static double betaLentz(double a, double b, double x) {
    constexpr int MaxIterations = 500;
    constexpr double MaxError = 1e-10;

    // Initialize coefficients.
    double sum = a + b;
    double q = a;
    double c = 1;
    double d = 1;
    updateD(d, -sum / (a + 1) * x);
    double frac = d;

    // In the continuous fraction expansion for the beta function,
    // odd and even coefficients are given by different formulas.
    // Therefore, each loop step below involves two updates.

    for (int m = 1; m <= MaxIterations; ++m) {
        double coeff;
        q += 2;
        coeff = m * (--b) / ((q - 1) * q);
        updateFraction(frac, c, d, coeff * x);
        coeff = (++a) * (++sum) / (q * (q + 1));
        double val = updateFraction(frac, c, d, -coeff * x) - 1;
        if (std::fabs(val) < MaxError)
            break;
    }
    return frac;
}

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

double zCDF(double z) {
    return zSF(-z);
}

// zSF(38) = 2.9e-316, zSF(39) = 0
double zSF(double z) {
    constexpr double SQRT2INV = 0.7071067811865475;  // 1 / sqrt(2)
    return 0.5 * std::erfc(z * SQRT2INV);
}

double tCDF(double t, double n) {
    double val = 0.5 * betainc(n / 2, 0.5, n / (n + t * t));
    return t < 0 ? val : 1 - val;
}

double tSF(double t, double n) {
    double val = 0.5 * betainc(n / 2, 0.5, n / (n + t * t));
    return t < 0 ? 1 - val : val;
}

double fCDF(double f, double n1, double n2) {
    return 1 - fSF(f, n1, n2);
}

double fSF(double f, double n1, double n2) {
    double x = n2 / (n2 + n1 * f);
    return betainc(n2 / 2, n1 / 2, x);
}

// Fast error function inverse using Winitzki's approximation.
// Expects 0 <= y <= 2.
// Could also compute p = std::log1p(-z * z) where z = y - 1,
// which is more accurate near y = 1 but overflows near 0 and 2.

double erfinv(double y) {
    constexpr double Pi = 3.141592653589793;
    constexpr double A = 7.1422296;
    constexpr double B = 2 * A / Pi;
    double p = std::log(y * (2 - y));
    double q = 0.5 * p + B;
    double r = std::sqrt(q * q - p * A) - q;
    return std::sqrt(2 * r);
}

double zCDFinv(double y) {
    auto val = erfinv(2 * y);
    return y >= 0.5 ? val : -val;
}

double t2z(double t, double n) {
    auto z = f2z(t * t, 1, n);
    return t >= 0 ? z : -z;
}

double f2z(double f, double n1, double n2) {
    return f > 0 ? erfinv(fSF(f, n1, n2)) : 0;
}

}  // namespace gxna
