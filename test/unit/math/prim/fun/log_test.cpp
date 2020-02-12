#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(MathFunctions, logInt) {
  EXPECT_FLOAT_EQ(std::log(3), stan::math::log(3));
  EXPECT_FLOAT_EQ(std::log(3.0), stan::math::log(3.0));
}

TEST(MathFunctions, log_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::log(b));
}
