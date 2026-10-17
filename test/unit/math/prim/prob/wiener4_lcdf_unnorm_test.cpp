#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>

namespace {

// log(Phi(-x) / phi(x)) against a 50-digit mpmath reference. The two terms
// grow like -x^2 / 2 while the ratio stays O(log x), so a difference of the
// two log densities loses digits without bound; the erfcx identity does
// not. The relative tolerance is 1e-14, about 45 ulp.
void expect_log_mill(double x, double expected) {
  const double got = stan::math::internal::logMill(x);
  EXPECT_NEAR(expected, got, 1e-14 * std::fabs(expected)) << "x=" << x;
}

}  // namespace

TEST(ProbWiener4LcdfUnnorm, logMill) {
  // negative and small arguments, where the difference form is exact
  expect_log_mill(-30, 450.91893853320467);
  expect_log_mill(-5, 13.41893824655306);
  expect_log_mill(-1, 1.2461847541812229);
  expect_log_mill(0, 0.22579135264472743);
  expect_log_mill(1, -0.42208311180459076);
  expect_log_mill(10, -2.3123466173077978);
  // large arguments, where the difference form loses 200 to 11 000 ulp
  expect_log_mill(30, -3.4023054231385244);
  expect_log_mill(50, -3.9124226062152285);
  expect_log_mill(100, -4.6052701610004159);
  expect_log_mill(300, -5.7037935854586871);
  expect_log_mill(1000, -6.9077562789796371);
}

TEST(ProbWiener4LcdfUnnorm, logMill_continuous_at_zero) {
  // the two evaluation forms meet at x = 0
  const double eps = 1e-12;
  const double left = stan::math::internal::logMill(-eps);
  const double right = stan::math::internal::logMill(eps);
  const double at = stan::math::internal::logMill(0.0);
  EXPECT_NEAR(at, left, 1e-11);
  EXPECT_NEAR(at, right, 1e-11);
}

TEST(ProbWiener4LcdfUnnorm, logMill_large_drift_finite) {
  // the public function at a point whose series reaches logMill with
  // arguments near 30: large drift, long time
  const double v = stan::math::wiener_lcdf_unnorm(100.0, 1.0, 0.0, 0.5, 3.0);
  EXPECT_TRUE(std::isfinite(v)) << v;
}
