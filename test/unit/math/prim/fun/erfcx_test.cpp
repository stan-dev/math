#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

// Reference values computed as expl(x * x) * erfcl(x) in long double
// (64-bit mantissa), which carries ~19 significant digits.
TEST(MathFunctions, erfcx) {
  using stan::math::erfcx;
  EXPECT_NEAR(108.94090438997797, erfcx(-2.0), 1.1e-12);
  EXPECT_NEAR(35.694955906025278, erfcx(-1.7), 3.6e-13);
  EXPECT_NEAR(5.0089800807622833, erfcx(-1.0), 5.1e-14);
  EXPECT_NEAR(1.952360489182557, erfcx(-0.5), 2.0e-14);
  EXPECT_NEAR(1.0, erfcx(0.0), 1e-14);
  EXPECT_NEAR(0.97782647768353936, erfcx(0.02), 1e-14);
  EXPECT_NEAR(0.6156903441929259, erfcx(0.5), 1e-14);
  EXPECT_NEAR(0.427583576155807, erfcx(1.0), 1e-14);
  EXPECT_NEAR(0.25539567631050575, erfcx(2.0), 1e-14);
  EXPECT_NEAR(0.15126529983237388, erfcx(3.6), 1e-14);
  EXPECT_NEAR(0.11070463773306863, erfcx(5.0), 1e-14);
}

// The whole point of the scaling: both factors of exp(x * x) * erfc(x) are
// out of range here, but the product is an ordinary number.
TEST(MathFunctions, erfcxUpperTail) {
  using stan::math::erfcx;
  EXPECT_EQ(0.0, std::erfc(30.0));  // the unscaled factor has underflowed
  EXPECT_NEAR(0.065925122499980351, erfcx(8.5), 1e-15);
  EXPECT_NEAR(0.046854221014893761, erfcx(12.0), 1e-15);
  EXPECT_NEAR(0.028174348741051319, erfcx(20.0), 1e-15);
  EXPECT_NEAR(0.021683584850562907, erfcx(26.0), 1e-15);
  EXPECT_NEAR(0.005641613782989433, erfcx(100.0), 1e-16);
  // collapses to the leading term 1 / (x * sqrt(pi)) once every correction
  // has fallen below the rounding of the result
  EXPECT_NEAR(stan::math::INV_SQRT_PI * 1e-10, erfcx(1e10), 1e-25);
}

TEST(MathFunctions, erfcxLowerTail) {
  using stan::math::erfcx;
  EXPECT_NEAR(144009798674.66104, erfcx(-5.0), 1.5e-3);
  EXPECT_NEAR(5.3762342836322712e+43, erfcx(-10.0), 5.4e+29);
  EXPECT_NEAR(7.6577249314905682e+293, erfcx(-26.0), 7.7e+279);
  // 2 * exp(x * x) leaves the binary64 range below about -26.63
  EXPECT_TRUE(std::isinf(erfcx(-27.0)));
  EXPECT_GT(erfcx(-27.0), 0.0);
}

// The internal crossover between exp(x * x) * erfc(x) and the Cody rational
// approximation must not be visible in the output.
TEST(MathFunctions, erfcxBranchContinuity) {
  using stan::math::erfcx;
  const double cut = 4.0;
  for (double delta : {1e-15, 1e-12, 1e-9, 1e-6}) {
    const double below = erfcx(cut - cut * delta);
    const double above = erfcx(cut + cut * delta);
    // the function is smooth and decreasing, so the two sides differ only by
    // the slope over 2 * cut * delta
    const double slope = 2.0 * cut * erfcx(cut) - stan::math::TWO_OVER_SQRT_PI;
    EXPECT_NEAR(below - above, -slope * 2.0 * cut * delta,
                1e-9 * delta + 1e-15);
  }
  EXPECT_NEAR(erfcx(std::nextafter(cut, 0.0)), erfcx(cut), 1e-16);
}

TEST(MathFunctions, erfcxEdgeCases) {
  using stan::math::erfcx;
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(1.0, erfcx(0.0));
  EXPECT_EQ(0.0, erfcx(inf));
  EXPECT_TRUE(std::isinf(erfcx(-inf)));
  EXPECT_GT(erfcx(-inf), 0.0);
  EXPECT_TRUE(std::isnan(erfcx(nan)));
}

TEST(MathFunctions, erfcxVectorized) {
  using stan::math::erfcx;
  std::vector<double> xs{-1.0, 0.0, 1.0, 20.0};
  std::vector<double> ys = erfcx(xs);
  ASSERT_EQ(4U, ys.size());
  for (size_t i = 0; i < xs.size(); ++i) {
    EXPECT_FLOAT_EQ(erfcx(xs[i]), ys[i]);
  }

  Eigen::VectorXd v(3);
  v << -2.0, 0.5, 12.0;
  Eigen::VectorXd w = erfcx(v);
  for (int i = 0; i < v.size(); ++i) {
    EXPECT_FLOAT_EQ(erfcx(v(i)), w(i));
  }
}

// erfcx is how the normal tail quantities stay finite; check the two
// identities that motivate the function.
TEST(MathFunctions, erfcxNormalTailIdentities) {
  using stan::math::erfcx;
  for (double x : {-40.0, -20.0, -7.5, -1.0}) {
    const double w = -x * stan::math::INV_SQRT_TWO;
    const double log_cdf = stan::math::LOG_HALF + std::log(erfcx(w)) - w * w;
    EXPECT_NEAR(log_cdf, stan::math::std_normal_lcdf(x),
                1e-10 * std::fabs(log_cdf));
    // inverse Mills ratio phi(x) / Phi(x)
    const double mills = stan::math::SQRT_TWO_OVER_SQRT_PI / erfcx(w);
    EXPECT_NEAR(
        mills,
        std::exp(stan::math::NEG_LOG_SQRT_TWO_PI - 0.5 * x * x - log_cdf),
        1e-9 * mills);
  }
}
