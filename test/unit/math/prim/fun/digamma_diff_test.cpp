#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, digamma_diff_special_cases) {
  using stan::math::digamma_diff;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();

  EXPECT_TRUE(std::isnan(digamma_diff(nan, 1.0)));
  EXPECT_TRUE(std::isnan(digamma_diff(1.0, nan)));
  EXPECT_EQ(digamma_diff(2.5, 0.0), 0.0);
  EXPECT_EQ(digamma_diff(1e18, 0), 0.0);
  EXPECT_EQ(digamma_diff(2.5, inf), inf);
  EXPECT_EQ(digamma_diff(inf, 2.5), 0.0);
  EXPECT_THROW(digamma_diff(0.0, 1.0), std::domain_error);
  EXPECT_THROW(digamma_diff(-1.5, 1.0), std::domain_error);
  EXPECT_THROW(digamma_diff(1.0, -0.5), std::domain_error);
}

TEST(MathFunctions, digamma_diff_integer_offset) {
  // psi(x + k) - psi(x) = sum_{j < k} 1 / (x + j) for integer k; the sum in
  // double has a rounding error of up to about k eps / 2
  using stan::math::digamma_diff;
  for (double x : {0.25, 1.0, 3.5, 9.5, 10.0, 47.0}) {
    double sum = 0;
    for (int k = 1; k <= 20; ++k) {
      sum += 1.0 / (x + k - 1);
      EXPECT_NEAR(digamma_diff(x, k), sum, 2e-15 * sum)
          << "x = " << x << ", k = " << k;
    }
  }
}

TEST(MathFunctions, digamma_diff_count_and_plain_paths) {
  // A count d from 0 to 8 is summed directly; for x < 10 and d >= 10 the
  // plain difference is used. d = 1 gives exactly 1 / x, as
  // digamma(x + 1) - digamma(x) does. References: mpmath at 80 digits.
  using stan::math::digamma_diff;
  for (double x : {1e-300, 0.1, 3.7, 9.99, 1e8, 1e300}) {
    EXPECT_EQ(digamma_diff(x, 1), 1.0 / x) << "x = " << x;
    EXPECT_EQ(digamma_diff(x, 1.0), 1.0 / x) << "x = " << x;
  }
  // count path
  EXPECT_NEAR(digamma_diff(2.5, 8), 1.599844393652443188, 1e-15 * 1.6);
  EXPECT_NEAR(digamma_diff(0x1.89374bc6a7efap-9, 7), 335.77886985018505786,
              1e-15 * 335.8);
  // plain difference, x < 10 and d >= 10
  EXPECT_NEAR(digamma_diff(0x1.ee45a1cac0831p+2, 10.0), 0.86831971186762265543,
              2e-15 * 0.87);
  EXPECT_NEAR(digamma_diff(0x1.7ae147ae147aep-2, 57.0), 6.8360822546693604884,
              2e-15 * 6.84);
  EXPECT_NEAR(digamma_diff(0x1.3fae147ae147bp+3, 1e6), 11.564819675088084883,
              2e-15 * 11.6);
  // d = 9 is neither a count of the first path nor large: the shift and the
  // asymptotic series
  EXPECT_NEAR(digamma_diff(6.25, 9.0), 0.9409809180725561924, 1e-15 * 0.94);
}

namespace digamma_diff_test_internal {
struct TestValue {
  double x;
  double d;
  double val;
};

// psi(x + d) - psi(x) computed with mpmath at 60 digits plus the digits lost
// to cancellation (log10(x) for large x, log10(1 / d) for small d); the
// arguments are written in hex so that they are exact. Points of large x
// are where the plain difference digamma(x + d) - digamma(x) is wrong: at
// x = 1e12, d = 57 it has a relative error of 4.6e-06, at x = 1e18 of 1.
std::vector<TestValue> testValues = {
    {0x1.0624dd2f1a9fcp-10, 0x1.0000000000000p-1, 9.9861698830979385e+2},
    {0x1.0000000000000p-1, 0x1.0000000000000p+0, 2.0000000000000000},
    {0x1.0000000000000p+0, 0x1.8000000000000p+1, 1.8333333333333333},
    {0x1.4000000000000p+1, 0x1.5798ee2308c3ap-27, 4.9035775491921462e-9},
    {0x1.3800000000000p+3, 0x1.0000000000000p-2, 2.6643054022145096e-2},
    {0x1.4000000000000p+3, 0x1.c800000000000p+5, 1.9454587802736456},
    {0x1.5000000000000p+3, 0x1.e848000000000p+19, 1.1512519523616630e+1},
    {0x1.2a00000000000p+5, 0x1.0000000000000p-1, 1.3512896711535966e-2},
    {0x1.f400000000000p+9, 0x1.d400000000000p+6, 1.1069890905638717e-1},
    {0x1.e848000000000p+19, 0x1.8000000000000p+1, 2.9999970000050000e-6},
    {0x1.7d78400000000p+26, 0x1.5798ee2308c3ap-27, 1.0000000050000000e-16},
    {0x1.d1a94a2000000p+39, 0x1.c800000000000p+5, 5.6999999998404000e-11},
    {0x1.b48eb57e00000p+44, 0x1.9000000000000p+8, 1.3333333333244667e-11},
    {0x1.c6bf526340000p+49, 0x1.e848000000000p+19, 9.9999999950000050e-10},
    {0x1.bc16d674ec800p+59, 0x1.0000000000000p+0, 1.0000000000000000e-18},
    {0x1.8232558201159p+59, 0x1.d400000000000p+6, 1.3453882098446929e-16},
    {0x1.0000000000000p-74, 0x1.0000000000000p+1, 1.8889465931478581e+22},
};
}  // namespace digamma_diff_test_internal

TEST(MathFunctions, digamma_diff_precomputed) {
  using digamma_diff_test_internal::TestValue;
  using digamma_diff_test_internal::testValues;
  using stan::math::digamma_diff;

  for (const TestValue& t : testValues) {
    EXPECT_NEAR(digamma_diff(t.x, t.d), t.val, 1e-15 * std::fabs(t.val))
        << "x = " << t.x << ", d = " << t.d;
  }
}
