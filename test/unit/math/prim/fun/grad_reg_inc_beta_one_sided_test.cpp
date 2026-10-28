#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>

namespace {

/**
 * The GRAD_RIB requests of grad_reg_inc_beta against the reference values
 * of grad_reg_inc_beta_test.cpp (mpmath, 80 digits). G1 and G2 carry only
 * their own derivative, with the other digamma unused. The comparison of
 * BOTH with the 9-argument form checks that the latter forwards to BOTH.
 */
void expect_one_sided(double a, double b, double z, double dda, double ddb,
                      double rtol) {
  using stan::math::digamma;
  using stan::math::GRAD_RIB;
  const double only_g1 = stan::math::grad_reg_inc_beta<GRAD_RIB::G1>(
      a, b, z, digamma(a), 0.0, digamma(a + b), 0.0);
  EXPECT_NEAR(dda, only_g1, rtol * std::fabs(dda))
      << "G1 a=" << a << " b=" << b << " z=" << z;
  const double only_g2 = stan::math::grad_reg_inc_beta<GRAD_RIB::G2>(
      a, b, z, 0.0, digamma(b), digamma(a + b), 0.0);
  EXPECT_NEAR(ddb, only_g2, rtol * std::fabs(ddb))
      << "G2 a=" << a << " b=" << b << " z=" << z;

  double g1 = 0;
  double g2 = 0;
  stan::math::grad_reg_inc_beta(g1, g2, a, b, z, digamma(a), digamma(b),
                                digamma(a + b), 0.0);
  const auto both = stan::math::grad_reg_inc_beta<GRAD_RIB::BOTH>(
      a, b, z, digamma(a), digamma(b), digamma(a + b), 0.0);
  EXPECT_EQ(g1, both.first) << "BOTH a=" << a << " b=" << b << " z=" << z;
  EXPECT_EQ(g2, both.second) << "BOTH a=" << a << " b=" << b << " z=" << z;
}

}  // namespace

TEST(grad_reg_inc_beta_one_sided, moderate) {
  expect_one_sided(1.0, 1.0, 0.4, -0.36651629274966202, 0.30649537425959442,
                   1e-12);
  expect_one_sided(2, 3, 0.25, -0.22805360232434637, 0.10692153005229139,
                   1e-12);
  expect_one_sided(0.5, 0.5, 0.5, -0.58312180806163756, 0.58312180806163756,
                   1e-12);
  expect_one_sided(20, 20, 0.5, -0.063740468830439859, 0.063740468830439859,
                   1e-12);
  expect_one_sided(1.5, 1.25, 0.6, -0.23806756196382869, 0.32279574755409757,
                   1e-12);
}

TEST(grad_reg_inc_beta_one_sided, underflow) {
  expect_one_sided(600, 600, 0.5, -0.01152047147304755, 0.01152047147304755,
                   1e-12);
  expect_one_sided(352, 590, 0.758, -3.3772717797123242e-138,
                   4.5390666074606649e-138, 1e-12);
}

TEST(grad_reg_inc_beta_one_sided, integer_b) {
  expect_one_sided(10, 1, 0.5, -0.00067690154351557159, 0.00226574342615196,
                   1e-12);
  expect_one_sided(1, 100, 0.01, -0.42882413065860014, 0.0036787479630394138,
                   1e-12);
  expect_one_sided(2500, 2, 0.999, -0.00020509953516671709, 0.24850317856904407,
                   1e-12);
}

TEST(grad_reg_inc_beta_one_sided, tails) {
  expect_one_sided(1.145, 6786, 1.8e-4, -0.38621356822721375,
                   5.8439109664263782e-5, 1e-10);
  expect_one_sided(446, 5, 0.76, -1.0644675765185877e-47,
                   1.2763067943072555e-46, 1e-12);
  expect_one_sided(15000, 1.25, 0.999, -6.5954322569556673e-10,
                   2.008497927784773e-6, 1e-10);
  expect_one_sided(15000, 12500, 0.5, -8.7210244281937301e-53,
                   9.5528279245811057e-53, 1e-10);
}
