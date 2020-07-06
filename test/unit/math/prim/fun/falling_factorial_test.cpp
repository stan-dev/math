#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, falling_factorial) {
  using stan::math::falling_factorial;

  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_FLOAT_EQ(24, falling_factorial(4.0, 3));
  EXPECT_FLOAT_EQ(120, falling_factorial(5.0, 4));
  EXPECT_THROW(falling_factorial(1, -4), std::domain_error);
  EXPECT_THROW(falling_factorial(nan, 1), std::domain_error);
  // see comments in test/unit/math/prim/fun/lgamma_test.cpp
  EXPECT_TRUE(std::isnormal(falling_factorial(1.0E30, 5)));
}

TEST(MathFunctions, falling_factorial_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::falling_factorial;
    return falling_factorial(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << -1.3, 0.7, 2.8;
  Eigen::VectorXi in2(3);
  in2 << 1, 3, 1;
  stan::test::binary_scalar_tester(f, in1, in2);
}
