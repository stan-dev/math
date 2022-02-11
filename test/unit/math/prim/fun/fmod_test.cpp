#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, fmodValues) {
  using stan::math::fmod;
  EXPECT_FLOAT_EQ(0.1, fmod(2.1, 2));
  EXPECT_FLOAT_EQ(1.9, fmod(10.0, 2.7));
  EXPECT_FLOAT_EQ(1.7, fmod(6, 4.3));
  EXPECT_FLOAT_EQ(0, fmod(6.0, 3.0));

  EXPECT_FLOAT_EQ(0, fmod(2, 1));
  EXPECT_FLOAT_EQ(0, fmod(2, 1.0));
  EXPECT_FLOAT_EQ(0, fmod(2.0, 1));
  EXPECT_FLOAT_EQ(0, fmod(2.0, 1.0));
}

TEST(MathFunctions, fmodNaN) {
  using stan::math::fmod;
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(fmod(1, nan)));
  EXPECT_TRUE(std::isnan(fmod(nan, 1)));
  EXPECT_TRUE(std::isnan(fmod(nan, nan)));
}

TEST(MathFunctions, fmodInf) {
  using stan::math::fmod;
  double inf = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(1, fmod(1, inf));
  EXPECT_TRUE(std::isnan(fmod(inf, 1)));
  EXPECT_TRUE(std::isnan(fmod(inf, inf)));
}

TEST(MathFunctions, fmod_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::fmod;
    return fmod(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
