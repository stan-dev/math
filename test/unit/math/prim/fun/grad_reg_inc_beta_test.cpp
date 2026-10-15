#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

namespace {

// Both shape derivatives of I_z(a, b) against a high-precision reference
// (mpmath, 80 digits, density-form quadrature checked against numerical
// differentiation at 110 digits).
void expect_grads(double a, double b, double z, double dda, double ddb,
                  double rtol) {
  using stan::math::digamma;
  using stan::math::lbeta;
  double g1 = 0;
  double g2 = 0;
  stan::math::grad_reg_inc_beta(g1, g2, a, b, z, digamma(a), digamma(b),
                                digamma(a + b), std::exp(lbeta(a, b)));
  EXPECT_TRUE(std::isfinite(g1)) << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_TRUE(std::isfinite(g2)) << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_NEAR(dda, g1, rtol * std::fabs(dda))
      << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_NEAR(ddb, g2, rtol * std::fabs(ddb))
      << "a=" << a << " b=" << b << " z=" << z;
}

}  // namespace

TEST(grad_reg_inc_beta, 1) {
  // at y = 1 the function is identically 1 and both derivatives are 0
  double alpha = 1.0;
  double beta = 1.0;
  double y = 1.0;
  double digamma_alpha = stan::math::digamma(alpha);
  double digamma_beta = stan::math::digamma(beta);
  double digamma_sum = stan::math::digamma(alpha + beta);
  double betafunc = std::exp(stan::math::lbeta(alpha, beta));

  double g1 = 0;
  double g2 = 0;
  stan::math::grad_reg_inc_beta(g1, g2, alpha, beta, y, digamma_alpha,
                                digamma_beta, digamma_sum, betafunc);
  EXPECT_FLOAT_EQ(0, g1);
  EXPECT_FLOAT_EQ(0, g2);
}

TEST(grad_reg_inc_beta, 2) {
  expect_grads(1.0, 1.0, 0.4, -0.36651629274966202, 0.30649537425959442, 1e-12);
}

TEST(grad_reg_inc_beta, moderate) {
  expect_grads(2, 3, 0.25, -0.22805360232434637, 0.10692153005229139, 1e-12);
  expect_grads(0.5, 0.5, 0.5, -0.58312180806163756, 0.58312180806163756, 1e-12);
  expect_grads(20, 20, 0.5, -0.063740468830439859, 0.063740468830439859, 1e-12);
  expect_grads(1.5, 1.25, 0.6, -0.23806756196382869, 0.32279574755409757,
               1e-12);
}

TEST(grad_reg_inc_beta, beta_ab_underflow) {
  // beta(a, b) underflows for a + b above about 1100
  expect_grads(600, 600, 0.5, -0.01152047147304755, 0.01152047147304755, 1e-12);
}

TEST(grad_reg_inc_beta, prefactor_underflow) {
  // z^a (1 - z)^b underflows while I is 1 and the derivatives are 1e-138
  expect_grads(352, 590, 0.758, -3.3772717797123242e-138,
               4.5390666074606649e-138, 1e-12);
}

TEST(grad_reg_inc_beta, integer_b) {
  // for integer b the power series of the value terminates but the
  // derivative series does not; the stop test must watch both
  expect_grads(10, 1, 0.5, -0.00067690154351557159, 0.00226574342615196, 1e-12);
  expect_grads(1, 100, 0.01, -0.42882413065860014, 0.0036787479630394138,
               1e-12);
  expect_grads(2500, 2, 0.999, -0.00020509953516671709, 0.24850317856904407,
               1e-12);
}

TEST(grad_reg_inc_beta, tails) {
  expect_grads(1.145, 6786, 1.8e-4, -0.38621356822721375, 5.8439109664263782e-5,
               1e-10);
  expect_grads(446, 5, 0.76, -1.0644675765185877e-47, 1.2763067943072555e-46,
               1e-12);
  expect_grads(15000, 1.25, 0.999, -6.5954322569556673e-10,
               2.008497927784773e-6, 1e-10);
  expect_grads(15000, 12500, 0.5, -8.7210244281937301e-53,
               9.5528279245811057e-53, 1e-10);
}
