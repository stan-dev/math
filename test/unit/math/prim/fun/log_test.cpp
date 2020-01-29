#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, logInt) {
  using stan::math::log;
  using std::log;
  EXPECT_FLOAT_EQ(std::log(3), log(3));
  EXPECT_FLOAT_EQ(std::log(3.0), log(3.0));
}

TEST(MathFunctions, log_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << 1.1, 1.2, 1.3, 1.4, 1.5;
  stan::math::multiply(a, stan::math::log(b));
}
