#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(MathFunctions, pow) {
  using stan::math::pow;
  using std::pow;

  EXPECT_FLOAT_EQ(pow(2.15, 1.71), 3.70228430936580996892756228068688);
  EXPECT_FLOAT_EQ(pow(7.62, -1.15), 0.0967724604023293179518727243992);
}

TEST(MathFunctions, pow_nan) {
  using stan::math::INFTY;
  using stan::math::NOT_A_NUMBER;
  using stan::math::pow;
  using std::pow;

  EXPECT_TRUE(std::isnan(pow(NOT_A_NUMBER, 2.16)));
  EXPECT_TRUE(std::isnan(pow(1.65, NOT_A_NUMBER)));
  EXPECT_TRUE(std::isnan(pow(-1.65, 2.16)));
  EXPECT_FALSE(std::isnan(pow(-1.65, 2.0)));

  EXPECT_TRUE(std::isinf(pow(INFTY, 2.16)));
  EXPECT_TRUE(std::isinf(pow(1.65, INFTY)));
  EXPECT_TRUE(std::isinf(pow(INFTY, INFTY)));
}

TEST(MathFunctions, pow_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using std::pow;
    using stan::math::pow;
    return pow(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.2, 3.1, 0.8;
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
