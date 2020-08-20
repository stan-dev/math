#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, hypotDouble) {
  using stan::math::hypot;

  EXPECT_FLOAT_EQ(5.0, hypot(3, 4));
  EXPECT_FLOAT_EQ(5.0, hypot(3, 4.0));
  EXPECT_FLOAT_EQ(5.0, hypot(3.0, 4));
  EXPECT_FLOAT_EQ(5.0, hypot(3.0, 4.0));

  EXPECT_FLOAT_EQ(0.0, hypot(0, 0));
  EXPECT_FLOAT_EQ(0.0, hypot(0, 0.0));
  EXPECT_FLOAT_EQ(0.0, hypot(0.0, 0));
  EXPECT_FLOAT_EQ(0.0, hypot(0.0, 0.0));
}

TEST(MathFunctions, hypotInf) {
  using stan::math::hypot;
  double inf = std::numeric_limits<double>::infinity();

  // promotes results to double
  EXPECT_FLOAT_EQ(inf, hypot(inf, 1));
  EXPECT_FLOAT_EQ(inf, hypot(1, inf));
  EXPECT_FLOAT_EQ(inf, hypot(inf, inf));

  EXPECT_FLOAT_EQ(inf, hypot(-inf, 1));
  EXPECT_FLOAT_EQ(inf, hypot(1, -inf));
  EXPECT_FLOAT_EQ(inf, hypot(-inf, -inf));
  EXPECT_FLOAT_EQ(inf, hypot(-inf, inf));
  EXPECT_FLOAT_EQ(inf, hypot(inf, -inf));
}

TEST(MathFunctions, hypotNaN) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::hypot(3.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::hypot(nan, 3.0)));

  EXPECT_TRUE(std::isnan(stan::math::hypot(nan, nan)));
}

TEST(MathFunctions, hypot_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::hypot;
    return hypot(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
