#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(MathFunctions, atan2) {
  using stan::math::atan2;

  EXPECT_FLOAT_EQ(atan2(2.15, 1.71), 0.8988979010770248000345);
  EXPECT_FLOAT_EQ(atan2(7.62, 10.15), 0.6439738474911284019668);
}

TEST(MathFunctions, atan2_nan) {
  using stan::math::atan2;
  using stan::math::INFTY;
  using stan::math::NOT_A_NUMBER;

  EXPECT_TRUE(std::isnan(atan2(NOT_A_NUMBER, 2.16)));
  EXPECT_TRUE(std::isnan(atan2(1.65, NOT_A_NUMBER)));

  EXPECT_FALSE(std::isnan(atan2(INFTY, 2.16)));
  EXPECT_FALSE(std::isnan(atan2(1.65, INFTY)));
}

TEST(MathFunctions, atan2_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::atan2(x1, x2);
  };

  Eigen::VectorXd in1 = Eigen::VectorXd::Random(6);
  Eigen::VectorXd in2 = Eigen::VectorXd::Random(6);

  stan::test::binary_scalar_tester(f, in1, in2);
}
