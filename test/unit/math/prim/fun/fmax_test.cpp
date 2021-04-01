#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, fmaxFinite) {
  using stan::math::fmax;
  EXPECT_FLOAT_EQ(1.0, fmax(1, 0));
  EXPECT_FLOAT_EQ(1.0, fmax(1.0, 0));
  EXPECT_FLOAT_EQ(1.0, fmax(1, 0.0));
  EXPECT_FLOAT_EQ(1.0, fmax(1.0, 0.0));

  EXPECT_FLOAT_EQ(1.0, fmax(0, 1));
  EXPECT_FLOAT_EQ(1.0, fmax(0, 1.0));
  EXPECT_FLOAT_EQ(1.0, fmax(0.0, 1));
  EXPECT_FLOAT_EQ(1.0, fmax(0.0, 1.0));
}

TEST(MathFunctions, fmaxNaN) {
  using stan::math::fmax;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FLOAT_EQ(1.0, fmax(1, nan));
  EXPECT_FLOAT_EQ(1.0, fmax(nan, 1));
  EXPECT_TRUE(std::isnan(stan::math::fmax(nan, nan)));
}

TEST(MathFunctions, fmaxInf) {
  using stan::math::fmax;
  double inf = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(inf, fmax(inf, 1));
  EXPECT_FLOAT_EQ(inf, fmax(1, inf));
  EXPECT_FLOAT_EQ(inf, fmax(inf, -inf));
  EXPECT_FLOAT_EQ(inf, fmax(-inf, inf));
}

TEST(MathFunctions, fmax_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::fmax;
    return fmax(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
