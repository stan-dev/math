#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/tgamma.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_gamma.hpp>

// converge
TEST(MathPrimScalFun, grad_reg_inc_gamma_1) {
  double alpha = 1.1;
  double z = 0.2;
  double g = stan::math::tgamma(alpha);
  double dig = stan::math::digamma(alpha);
  EXPECT_NEAR(0.31416364892410884,
              stan::math::grad_reg_inc_gamma(alpha, z, g, dig, 1e-10), 1e-8);
}
// converge
TEST(MathPrimScalFun, grad_reg_inc_gamma_2) {
  double alpha = 1.1;
  double z = 2.0;
  double g = stan::math::tgamma(alpha);
  double dig = stan::math::digamma(alpha);
  EXPECT_NEAR(0.2350546737889920,
              stan::math::grad_reg_inc_gamma(alpha, z, g, dig, 1e-10), 1e-8);
}
// converge
TEST(MathPrimScalFun, grad_reg_inc_gamma_3) {
  double alpha = 2.5;
  double z = 1.3;
  double g = stan::math::tgamma(alpha);
  double dig = stan::math::digamma(alpha);
  EXPECT_NEAR(0.22962689833555939,
              stan::math::grad_reg_inc_gamma(alpha, z, g, dig, 1e-10), 1e-8);
}
// converge
TEST(MathPrimScalFun, grad_reg_inc_gamma_4) {
  double alpha = 2.5;
  double z = 30;
  double g = stan::math::tgamma(alpha);
  double dig = stan::math::digamma(alpha);
  EXPECT_NEAR(3.3205e-11,
              stan::math::grad_reg_inc_gamma(alpha, z, g, dig, 1e-10), 1e-8);
}
// converge
TEST(MathPrimScalFun, grad_reg_inc_gamma_5) {
  double alpha = 9;
  double z = 10;
  double g = stan::math::tgamma(alpha);
  double dig = stan::math::digamma(alpha);
  EXPECT_NEAR(0.120855166827777,
              stan::math::grad_reg_inc_gamma(alpha, z, g, dig, 1e-10), 4e-7);
}
// converge
TEST(MathPrimScalFun, grad_reg_inc_gamma_6) {
  double alpha = 10.0;
  double z = 9.0;
  double g = stan::math::tgamma(alpha);
  double dig = stan::math::digamma(alpha);
  EXPECT_NEAR(0.1270365119242684,
              stan::math::grad_reg_inc_gamma(alpha, z, g, dig, 1e-12), 1e-8);
}
