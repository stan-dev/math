#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, inc_beta_inv_inv) {
  using stan::math::inc_beta_inv;

  EXPECT_FLOAT_EQ(0.146446609407, inc_beta_inv(0.5, 0.5, 0.25));
  EXPECT_FLOAT_EQ(0.5, inc_beta_inv(0.5, 0.5, 0.5));
  EXPECT_FLOAT_EQ(0.853553390593, inc_beta_inv(0.5, 0.5, 0.75));
  EXPECT_FLOAT_EQ(1.0, inc_beta_inv(0.5, 0.5, 1.0));

  EXPECT_FLOAT_EQ(0.0, inc_beta_inv(0.1, 1.5, 0.0));
  EXPECT_FLOAT_EQ(5.33624995327e-7, inc_beta_inv(0.1, 1.5, 0.25));
  EXPECT_FLOAT_EQ(0.000546567646391, inc_beta_inv(0.1, 1.5, 0.5));
  EXPECT_FLOAT_EQ(0.0319736161201, inc_beta_inv(0.1, 1.5, 0.75));
  EXPECT_FLOAT_EQ(1.0, inc_beta_inv(0.1, 1.5, 1.0));

  EXPECT_FLOAT_EQ(0.804032164227, inc_beta_inv(0.6, 0.3, 0.5));
}

TEST(MathFunctions, inc_beta_inv_a_boundary) {
  double b = 0.5;
  double p = 0.5;

  const double inf = std::numeric_limits<double>::infinity();

  EXPECT_THROW(stan::math::inc_beta_inv(0.0, b, p), std::domain_error);
  EXPECT_NO_THROW(stan::math::inc_beta_inv(inf, b, p));
  EXPECT_THROW(stan::math::inc_beta_inv(-0.01, b, p), std::domain_error);
}

TEST(MathFunctions, inc_beta_inv_b_boundary) {
  double a = 0.5;
  double p = 0.5;

  const double inf = std::numeric_limits<double>::infinity();

  EXPECT_THROW(stan::math::inc_beta_inv(a, 0.0, p), std::domain_error);
  EXPECT_NO_THROW(stan::math::inc_beta_inv(inf, 1.0, p));
  EXPECT_THROW(stan::math::inc_beta_inv(a, -0.01, p), std::domain_error);
}

TEST(MathFunctions, inc_beta_inv_x_boundary) {
  double a = 0.5;
  double b = 0.5;

  EXPECT_NO_THROW(stan::math::inc_beta_inv(a, b, 0.0));
  EXPECT_NO_THROW(stan::math::inc_beta_inv(a, b, 1.0));
  EXPECT_THROW(stan::math::inc_beta_inv(a, b, -0.01), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta_inv(a, b, 1.01), std::domain_error);
}

TEST(MathFunctions, inc_beta_inv_a_b_boundary) {
  double p = 0.5;

  EXPECT_THROW(stan::math::inc_beta_inv(0.0, 0.0, p), std::domain_error);
}

TEST(MathFunctions, inc_beta_inv_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::inc_beta_inv(0.5, 0.0, nan), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta_inv(0.5, nan, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta_inv(0.5, nan, nan), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta_inv(nan, 0.0, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta_inv(nan, 0.0, nan), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta_inv(nan, nan, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta_inv(nan, nan, nan), std::domain_error);
}
