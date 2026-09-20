#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>

namespace {

// d/da of I_z(a, b) against a high-precision reference (mpmath, 80 digits,
// density-form quadrature checked against numerical differentiation).
void expect_dda(double a, double b, double z, double expected, double rtol) {
  using stan::math::digamma;
  using stan::math::inc_beta_dda;
  const double got = inc_beta_dda(a, b, z, digamma(a), digamma(a + b));
  EXPECT_TRUE(std::isfinite(got)) << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_NEAR(expected, got, rtol * std::fabs(expected))
      << "a=" << a << " b=" << b << " z=" << z;
}

}  // namespace

TEST(MathFunctions, inc_beta_dda) {
  using stan::math::digamma;
  using stan::math::inc_beta_dda;

  double small_a = 1.5;
  double large_a = 15000;

  double small_b = 1.25;
  double large_b = 12500;

  double small_z = 0.001;
  double mid_z = 0.6;
  double large_z = 0.999;

  expect_dda(small_a, small_b, small_z, -0.00028665636570426377, 1e-10);
  expect_dda(small_a, small_b, mid_z, -0.23806756196382869, 1e-10);
  expect_dda(small_a, small_b, large_z, -0.00022264492614285471, 1e-10);

  // I underflows to 0 in double; the derivative is below 1e-3000
  EXPECT_FLOAT_EQ(0.0, inc_beta_dda(large_a, small_b, small_z, digamma(large_a),
                                    digamma(large_a + small_b)));
  EXPECT_FLOAT_EQ(0.0, inc_beta_dda(large_a, small_b, mid_z, digamma(large_a),
                                    digamma(large_a + small_b)));
  expect_dda(large_a, small_b, large_z, -6.5954322569556673e-10, 1e-10);

  expect_dda(small_a, large_b, small_z, -3.937563719959171e-5, 1e-10);
  EXPECT_FLOAT_EQ(0.0, inc_beta_dda(small_a, large_b, mid_z, digamma(small_a),
                                    digamma(small_a + large_b)));
  EXPECT_FLOAT_EQ(0.0, inc_beta_dda(small_a, large_b, large_z, digamma(small_a),
                                    digamma(small_a + large_b)));

  EXPECT_FLOAT_EQ(0.0, inc_beta_dda(large_a, large_b, small_z, digamma(large_a),
                                    digamma(large_a + large_b)));
  expect_dda(large_a, large_b, mid_z, -1.755811372560583e-76, 1e-10);
  EXPECT_FLOAT_EQ(0.0, inc_beta_dda(large_a, large_b, large_z, digamma(large_a),
                                    digamma(large_a + large_b)));
}

TEST(MathFunctions, inc_beta_dda_regressions) {
  // b >> a with small z
  expect_dda(1.145, 6786, 1.8e-4, -0.38621356822721375, 1e-10);
  // deep lower tail, I of order 1e-47
  expect_dda(446, 5, 0.76, -1.0644675765185877e-47, 1e-10);
  // large a with z near 1
  expect_dda(2500, 2, 0.999, -0.00020509953516671709, 1e-10);
  expect_dda(1, 100, 0.01, -0.42882413065860014, 1e-10);
  // integer b: the value series terminates but the derivative does not
  expect_dda(2, 3, 0.25, -0.22805360232434637, 1e-10);
  expect_dda(10, 1, 0.5, -0.00067690154351557159, 1e-10);
  // a + b above the beta(a, b) underflow limit
  expect_dda(600, 600, 0.5, -0.01152047147304755, 1e-10);
  expect_dda(0.5, 0.5, 0.5, -0.58312180806163756, 1e-10);
  expect_dda(20, 20, 0.5, -0.063740468830439859, 1e-10);
}
