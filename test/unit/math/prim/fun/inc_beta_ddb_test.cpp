#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>

namespace {

// d/db of I_z(a, b) against a high-precision reference (mpmath, 80 digits,
// density-form quadrature checked against numerical differentiation).
// The fourth argument of inc_beta_ddb is digamma(b).
void expect_ddb(double a, double b, double z, double expected, double rtol) {
  using stan::math::digamma;
  using stan::math::inc_beta_ddb;
  const double got = inc_beta_ddb(a, b, z, digamma(b), digamma(a + b));
  EXPECT_TRUE(std::isfinite(got)) << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_NEAR(expected, got, rtol * std::fabs(expected))
      << "a=" << a << " b=" << b << " z=" << z;
}

}  // namespace

TEST(MathFunctions, inc_beta_ddb) {
  using stan::math::digamma;
  using stan::math::inc_beta_ddb;

  double small_a = 1.5;
  double large_a = 15000;

  double small_b = 1.25;
  double large_b = 12500;

  double small_z = 0.001;
  double mid_z = 0.5;
  double large_z = 0.999;

  expect_ddb(small_a, small_b, small_z, 4.4135732827365698e-5, 1e-10);
  expect_ddb(small_a, small_b, mid_z, 0.29301795311347039, 1e-10);
  expect_ddb(small_a, small_b, large_z, 0.0018969609780283825, 1e-10);

  // I underflows to 0 in double; the derivative is below 1e-3000
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddb(large_a, small_b, small_z, digamma(small_b),
                                    digamma(large_a + small_b)));
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddb(large_a, small_b, mid_z, digamma(small_b),
                                    digamma(large_a + small_b)));
  expect_ddb(large_a, small_b, large_z, 2.008497927784773e-6, 1e-10);

  expect_ddb(small_a, large_b, small_z, 1.4782043391986094e-8, 1e-10);
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddb(small_a, large_b, mid_z, digamma(large_b),
                                    digamma(small_a + large_b)));
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddb(small_a, large_b, large_z, digamma(large_b),
                                    digamma(small_a + large_b)));

  EXPECT_FLOAT_EQ(0.0, inc_beta_ddb(large_a, large_b, small_z, digamma(large_b),
                                    digamma(large_a + large_b)));
  expect_ddb(large_a, large_b, mid_z, 9.5528279245811057e-53, 1e-10);
  EXPECT_FLOAT_EQ(0.0, inc_beta_ddb(large_a, large_b, large_z, digamma(large_b),
                                    digamma(large_a + large_b)));
}

TEST(MathFunctions, inc_beta_ddb_regressions) {
  expect_ddb(1.145, 6786, 1.8e-4, 5.8439109664263782e-5, 1e-10);
  expect_ddb(446, 5, 0.76, 1.2763067943072555e-46, 1e-10);
  expect_ddb(2500, 2, 0.999, 0.24850317856904407, 1e-10);
  expect_ddb(1, 100, 0.01, 0.0036787479630394138, 1e-10);
  // integer b: the value series terminates but the derivative does not
  expect_ddb(2, 3, 0.25, 0.10692153005229139, 1e-10);
  expect_ddb(10, 1, 0.5, 0.00226574342615196, 1e-10);
  expect_ddb(600, 600, 0.5, 0.01152047147304755, 1e-10);
  expect_ddb(0.5, 0.5, 0.5, 0.58312180806163756, 1e-10);
  expect_ddb(20, 20, 0.5, 0.063740468830439859, 1e-10);
}
