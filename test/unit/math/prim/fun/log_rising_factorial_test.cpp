#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, log_rising_factorial) {
  using stan::math::log_rising_factorial;

  EXPECT_FLOAT_EQ(std::log(120.0), log_rising_factorial(4.0, 3));
  EXPECT_FLOAT_EQ(std::log(360.0), log_rising_factorial(3.0, 4));
  EXPECT_THROW(log_rising_factorial(-1, 4), std::domain_error);
}

TEST(MathFunctions, log_rising_factorial_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log_rising_factorial(nan, 3)));
}

TEST(MathFunctions, log_rising_factorial_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_rising_factorial;
    return log_rising_factorial(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << 7.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
