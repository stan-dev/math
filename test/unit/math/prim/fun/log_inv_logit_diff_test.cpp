#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(MathFunctions, log_inv_logit_diff_inf) {
  using stan::math::INFTY;
  using stan::math::log_inv_logit_diff;

  EXPECT_FLOAT_EQ(-0.0916494, log_inv_logit_diff(INFTY, -2.34361));
  EXPECT_TRUE(std::isnan(log_inv_logit_diff(-INFTY, -2.34361)));
}

TEST(MathFunctions, log_inv_logit_diff) {
  using stan::math::log_inv_logit_diff;

  EXPECT_FLOAT_EQ(-3.019359355, log_inv_logit_diff(2.15, 1.71));
  EXPECT_FLOAT_EQ(-7.703540544, log_inv_logit_diff(-7.62, -10.15));
}

TEST(MathFunctions, log_inv_logit_diff_nan) {
  using stan::math::log_inv_logit_diff;
  using stan::math::NOT_A_NUMBER;

  EXPECT_TRUE(std::isnan(log_inv_logit_diff(NOT_A_NUMBER, 2.16)));
}

TEST(MathFunctions, log_inv_logit_diff_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_inv_logit_diff;
    return log_inv_logit_diff(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 4.2;
  Eigen::VectorXd in2(3);
  in2 << 0.3, 0.7, -2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
