#include "Distribution.h"

#include <cmath>
#include "gtest/gtest.h"

namespace gxna {

TEST(Distribution, BetaInc) {
    // Special values
    double eps = 1e-8;
    for (double x = 0; x <= 1; x += 0.01) {
        EXPECT_NEAR(betainc(1, 1, x), x, eps);
        EXPECT_NEAR(betainc(2, 1, x), x * x, eps);
        EXPECT_NEAR(betainc(4, 1, x), x * x * x * x, eps);
        EXPECT_NEAR(betainc(1, 2, x), x * (2 - x), eps);
    }
}

TEST(Distribution, zCDF) {
    // Asymptotics
    EXPECT_DOUBLE_EQ(zCDF(-40), 0);
    EXPECT_NEAR(zCDF(-20), 0, 1e-80);
    EXPECT_NEAR(zSF(20), 0, 1e-80);
    EXPECT_DOUBLE_EQ(zSF(40), 0);

    // Special values
    double eps = 1e-9;
    EXPECT_DOUBLE_EQ(zCDF(0), 0.5);
    EXPECT_NEAR(zCDF(1), 0.841344746, eps);
    EXPECT_NEAR(zCDF(2), 0.977249868, eps);
    EXPECT_NEAR(zCDF(3), 0.998650102, eps);
}

TEST(Distribution, tCDF) {
    // Asymptotics
    for (double n = 0.2; n <= 10; ++n) {
        EXPECT_NEAR(tCDF(-1e20, n), 0, 1e-4);
        EXPECT_EQ(tCDF(0, n), 0.5);
    }

    // Special values
    constexpr double Pi = 3.141592653589793;
    for (double t = -10; t <= 10; t += 0.2) {
        EXPECT_NEAR(tCDF(t, 1), 0.5 + std::atan(t) / Pi, 1e-9);  // t(1) is Cauchy
        EXPECT_NEAR(tCDF(t, 2), 0.5 + t / (2 * std::sqrt(2 + t * t)), 1e-9);
        EXPECT_NEAR(tCDF(t, 1e6), zCDF(t), 1e-6);  // t(inf) is Normal
    }
}

// F-distribution tail asymptotics
double fSFApprox(double f, double n1, double n2) {
    n1 /= 2;
    n2 /= 2;
    double beta = std::tgamma(n1) * std::tgamma(n2) / std::tgamma(n1 + n2);
    double coeff = 1 / (beta * std::pow(n1 / n2, n2) * n2);
    return std::pow(f, -n2) * coeff;
}

TEST(Distribution, fCDF) {
    double eps = 1e-9;

    for (double n1 = 0.5; n1 <= 10; n1 += 0.2)
        for (double n2 = 0.5; n2 <= 10; n2 += 0.2) {
            // Zero to zero
            EXPECT_EQ(fCDF(0, n1, n2), 0);

            // Asymptotics
            for (auto f : { 1e4, 1e5, 1e6, 1e7 }) {
                double val = fSF(f, n1, n2) / fSFApprox(f, n1, n2);
                EXPECT_NEAR(val, 1, 100 / f) << f << ' ' << n1 << ' ' << n2;
            }

            // If X ~ F(n1, n2) then 1 / X ~ F(n2, n1)
            for (auto x : { 0.3, 0.4, 0.5, 2.0, 3.0, 4.0 })
                EXPECT_NEAR(fCDF(x, n1, n2), 1 - fCDF(1 / x, n2, n1), eps);
        }

    // Symmetry
    for (double n = 0.1; n <= 10; n += 0.2)
        EXPECT_NEAR(fCDF(1, n, n), 0.5, eps);

    // Special values
    EXPECT_NEAR(fCDF(3.0, 1, 1), 2 / 3., eps);
    EXPECT_NEAR(fCDF(3.0, 2, 2), 3 / 4., eps);
    EXPECT_NEAR(fCDF(5.0, 10, 20), .998903410, eps);
}

TEST(Distribution, T2Z) {
    // T(Inf) is Normal
    for (double t = -38; t <= 38; t += 1)
        EXPECT_NEAR(t, t2z(t, 1e5), 0.005 * std::fabs(t));
}

TEST(Distribution, zCDFInv) {
    // Symmetry
    EXPECT_EQ(zCDFinv(0.5), 0);

    // Special values (note absolute error is fairly high)
    EXPECT_NEAR(zCDFinv(1e-316), -38.0278, 0.04);
    EXPECT_NEAR(zCDFinv(1e-80), -18.9916, 0.04);
    EXPECT_NEAR(zCDFinv(1e-40), -13.3109, 0.04);
    EXPECT_NEAR(zCDFinv(1e-20), -9.2623, 0.015);
    EXPECT_NEAR(zCDFinv(1e-10), -6.3613, 0.01);
    EXPECT_NEAR(zCDFinv(0.01), -2.3263, 0.01);
    EXPECT_NEAR(zCDFinv(0.025), -1.9599, 0.01);
    EXPECT_NEAR(zCDFinv(0.8414), 1.0002, 0.0002);
}

}  // namespace gxna
