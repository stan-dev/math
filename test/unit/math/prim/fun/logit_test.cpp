#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, logit) {
  using stan::math::logit;
  EXPECT_FLOAT_EQ(0.0, logit(0.5));
  EXPECT_FLOAT_EQ(5.0, logit(1.0 / (1.0 + exp(-5.0))));
}

TEST(MathFunctions, logit_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::logit(nan)));
}

TEST(MathFunctions, logit_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::logit(b));
}
