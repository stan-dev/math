#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(MathFunctions, pow) {
  EXPECT_FLOAT_EQ(stan::math::pow(2, 7), std::pow(2, 7));
  EXPECT_FLOAT_EQ(stan::math::pow(2.2, 7), std::pow(2.2, 7));
  EXPECT_FLOAT_EQ(stan::math::pow(2.15, 1.71), std::pow(2.15, 1.71));
  EXPECT_FLOAT_EQ(stan::math::pow(7.62, -1.15), std::pow(7.62, -1.15));
}

TEST(MathFunctions, pow_nan) {
  using stan::math::INFTY;
  using stan::math::NOT_A_NUMBER;
  using stan::math::pow;

  EXPECT_TRUE(std::isnan(pow(NOT_A_NUMBER, 2.16)));
  EXPECT_TRUE(std::isnan(pow(1.65, NOT_A_NUMBER)));
  EXPECT_TRUE(std::isnan(pow(-1.65, 2.16)));
  EXPECT_FALSE(std::isnan(pow(-1.65, 2.0)));

  EXPECT_TRUE(std::isinf(pow(INFTY, 2.16)));
  EXPECT_TRUE(std::isinf(pow(1.65, INFTY)));
  EXPECT_TRUE(std::isinf(pow(INFTY, INFTY)));
}

TEST(MathFunctions, pow_vec) {
  auto f
      = [](const auto& x1, const auto& x2) { return stan::math::pow(x1, x2); };

  Eigen::VectorXd in1(3);
  in1 << 1.2, 3.1, 0.8;
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
