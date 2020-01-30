#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log_inv_logit) {
  using stan::math::inv_logit;
  using stan::math::log_inv_logit;
  using std::log;

  EXPECT_FLOAT_EQ(log(inv_logit(-7.2)), log_inv_logit(-7.2));
  EXPECT_FLOAT_EQ(log(inv_logit(0.0)), log_inv_logit(0.0));
  EXPECT_FLOAT_EQ(log(inv_logit(1.9)), log_inv_logit(1.9));
}

TEST(MathFunctions, log_inv_logit_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log_inv_logit(nan)));
}

TEST(MathFunctions, log_inv_logit_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::log_inv_logit(b));
}
