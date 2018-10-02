#include <stan/math/prim/scal/fun/log_inv_logit_diff.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, log_inv_logit_diff) {
  using stan::math::log_inv_logit_diff;

  EXPECT_FLOAT_EQ(-3.019359355, log_inv_logit_diff(2.15, 1.71));
  EXPECT_FLOAT_EQ(-7.703540544, log_inv_logit_diff(-7.62, -10.15));
}

TEST(MathFunctions, log_inv_logit_diff_nan) {
  using stan::math::NOT_A_NUMBER;
  using stan::math::log_inv_logit_diff;

  EXPECT_PRED1(boost::math::isnan<double>,
               log_inv_logit_diff(NOT_A_NUMBER, 2.16));
}
