#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, sqrtInt) {
  EXPECT_FLOAT_EQ(std::sqrt(3.0), stan::math::sqrt(3));
  EXPECT_TRUE(stan::math::is_nan(stan::math::sqrt(-2)));
}

TEST(MathFunctions, sqrt_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::sqrt(b));
}
