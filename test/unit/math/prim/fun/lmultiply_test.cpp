#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, lmultiply) {
  double a = 2.0;
  double b = 3.0;
  EXPECT_FLOAT_EQ(a * log(b), stan::math::lmultiply(a, b));

  a = 0.0;
  b = 0.0;
  EXPECT_FLOAT_EQ(0.0, stan::math::lmultiply(a, b))
      << "when a and b are both 0, the result should be 0";

  a = 1.0;
  b = -1.0;
  EXPECT_TRUE(std::isnan(stan::math::lmultiply(a, b)))
      << "log(b) with b < 0 should result in NaN";
}

TEST(MathFunctions, lmultiply_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::lmultiply(2.0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::lmultiply(nan, 3.0)));
  EXPECT_TRUE(std::isnan(stan::math::lmultiply(nan, nan)));
}

TEST(MathFunctions, lmultiply_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::lmultiply;
    return lmultiply(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << 7.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
