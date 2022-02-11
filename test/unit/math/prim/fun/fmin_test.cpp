#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, fminFinite) {
  using stan::math::fmin;
  EXPECT_FLOAT_EQ(1.0, fmin(1, 2));
  EXPECT_FLOAT_EQ(1.0, fmin(1.0, 2));
  EXPECT_FLOAT_EQ(1.0, fmin(1, 2.0));
  EXPECT_FLOAT_EQ(1.0, fmin(1.0, 2.0));

  EXPECT_FLOAT_EQ(1.0, fmin(2, 1));
  EXPECT_FLOAT_EQ(1.0, fmin(2, 1.0));
  EXPECT_FLOAT_EQ(1.0, fmin(2.0, 1));
  EXPECT_FLOAT_EQ(1.0, fmin(2.0, 1.0));
}

TEST(MathFunctions, fminNaN) {
  using stan::math::fmin;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FLOAT_EQ(1.0, fmin(1, nan));
  EXPECT_FLOAT_EQ(1.0, fmin(nan, 1));
  EXPECT_TRUE(std::isnan(stan::math::fmin(nan, nan)));
}

TEST(MathFunctions, fminInf) {
  using stan::math::fmin;
  double inf = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(1, fmin(inf, 1));
  EXPECT_FLOAT_EQ(1, fmin(1, inf));
  EXPECT_FLOAT_EQ(-inf, fmin(inf, -inf));
  EXPECT_FLOAT_EQ(-inf, fmin(-inf, inf));
}

TEST(MathFunctions, fmin_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::fmin;
    return fmin(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
