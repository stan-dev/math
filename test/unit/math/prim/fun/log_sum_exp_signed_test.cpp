#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <tuple>
#include <vector>

TEST(MathFunctions, log_sum_exp_signed) {
  using stan::math::exp;
  using stan::math::log_diff_exp;
  using stan::math::log_sum_exp;
  using stan::math::log_sum_exp_signed;
  using stan::math::sign;

  double a = 2.5;
  double b = 76.2;
  double exp_a = exp(a);
  double exp_b = exp(b);

  double val;
  int val_sign;

  std::forward_as_tuple(val, val_sign) = log_sum_exp_signed(a, 1, b, 1);

  EXPECT_FLOAT_EQ(exp_a + exp_b, val_sign * exp(val));

  std::forward_as_tuple(val, val_sign) = log_sum_exp_signed(a, 1, b, -1);

  EXPECT_FLOAT_EQ(exp_a - exp_b, val_sign * exp(val));

  std::forward_as_tuple(val, val_sign) = log_sum_exp_signed(a, -1, b, 1);

  EXPECT_FLOAT_EQ(-exp_a + exp_b, val_sign * exp(val));

  std::forward_as_tuple(val, val_sign) = log_sum_exp_signed(a, -1, b, -1);

  EXPECT_FLOAT_EQ(-exp_a - exp_b, val_sign * exp(val));
}
