#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, inv_logit) {
  using stan::math::inv_logit;
  EXPECT_FLOAT_EQ(0.5, inv_logit(0.0));
  EXPECT_FLOAT_EQ(1.0 / (1.0 + exp(-5.0)), inv_logit(5.0));
}

TEST(MathFunctions, inv_logit_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::inv_logit(nan)));
}

TEST(MathFunctions, inv_logit_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::inv_logit(b));
}
