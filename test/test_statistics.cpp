#include "Statistics.h"

#include <cmath>
#include <limits>
#include <vector>
#include "gtest/gtest.h"

namespace gxna {

// Data sets for testing
class DataSetTest : public testing::Test {
 protected:
    DataSetTest()
        : fa(a), fb(b), fc(c)
    {}

    const std::vector<double> a { 1, 3, 2, 4, 5, 2.5, 0, 1, 0.8 };
    const std::vector<double> b { 5, 3, 8, 8.1, 10, 7, 7, 9, 0.2, 4, 4, 4.5 };
    const std::vector<double> c { 3, 2, 7, 3.1, 8, 9, 0, 0 };

    FastDataSet fa;
    FastDataSet fb;
    FastDataSet fc;
};

TEST_F(DataSetTest, Basic) {
    FastDataSet fz;
    fz += fa;
    fz += fb;
    fz += fc;

    EXPECT_TRUE(fz.size() == fa.size() + fb.size() + fc.size());

    EXPECT_DOUBLE_EQ(fz.mean(), 4.179310344827587);
    EXPECT_DOUBLE_EQ(fz.sd(), 3.0392126147209884);
}

TEST_F(DataSetTest, TF) {
    // Special values and symmetry
    double t = 3.7245441941478252;
    EXPECT_DOUBLE_EQ(tstat(fa, fb), t);
    EXPECT_DOUBLE_EQ(tstat(fb, fa), -t);

    double tequal = 3.454821662116653;
    EXPECT_DOUBLE_EQ(tstatEqual(fa, fb), tequal);
    EXPECT_DOUBLE_EQ(tstatEqual(fb, fa), -tequal);

    EXPECT_DOUBLE_EQ(fstat({fa, fb}), 11.935792717030486);
    EXPECT_DOUBLE_EQ(fstat({fa, fb, fc}), 4.5688416461171313);

    constexpr auto Inf = std::numeric_limits<double>::infinity();

    // t stat corner cases
    EXPECT_EQ(tstat(fa, fa), 0);
    EXPECT_EQ(tstat(fb, fb), 0);
    FastDataSet ones {{1, 1, 1, 1}};
    FastDataSet twos {{2, 2, 2, 2, 2, 2}};
    EXPECT_EQ(tstat(ones, ones), 0);
    EXPECT_EQ(tstat(twos, twos), 0);
    EXPECT_EQ(tstat(ones, twos), Inf);
    EXPECT_EQ(tstat(twos, ones), -Inf);

    // F stat corner cases
    EXPECT_EQ(fstat({fa, fa}), 0);
    EXPECT_EQ(fstat({fb, fb}), 0);
    EXPECT_EQ(fstat({ones, ones}), 0);
    EXPECT_EQ(fstat({twos, twos}), 0);
    EXPECT_EQ(fstat({ones, twos}), Inf);
    EXPECT_EQ(fstat({twos, ones}), Inf);

    // Singletons
    EXPECT_EQ(fstat({ FastDataSet({1}), FastDataSet({2}), FastDataSet({5}) }), 0);
}

TEST(Gamma, Basic) {
    constexpr double Pi = 3.141592653589793;
    constexpr double Gamma = 0.577215664901532;
    double eps = 1e-9;

    // Spot checks
    EXPECT_NEAR(harmonic(20), 3.5977396571436, eps);
    EXPECT_NEAR(harmonic2(20), 1.59616324391302, eps);

    // Special values
    EXPECT_NEAR(digamma(0.5), -std::log(4) - Gamma, eps);
    EXPECT_NEAR(digamma(1), -Gamma, eps);
    EXPECT_NEAR(digamma(2), 1 - Gamma, eps);

    // Spot checks
    EXPECT_NEAR(digamma(10), 2.25175258906, eps);
    EXPECT_NEAR(digamma(100), 4.60016185273, eps);

    // Special values
    EXPECT_NEAR(trigamma(0.5), Pi * Pi / 2, eps);
    EXPECT_NEAR(trigamma(1), Pi * Pi / 6, eps);

    // Spot checks
    EXPECT_NEAR(trigamma(10), 0.105166335, eps);
    EXPECT_NEAR(trigamma(100), 0.010050167, eps);
}

TEST(Gamma, Inverse) {
    for (double x = 1e-10; x <= 1e10; x *= 2) {
        double y = trigamma(x);
        EXPECT_NEAR(x, trigammainv(y), 1e-6);
    }
}

}  // namespace gxna
