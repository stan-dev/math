#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, squared_distance) {
  double x1 = 1;
  double x2 = 4;

  EXPECT_FLOAT_EQ(9, stan::math::squared_distance(x1, x2));
  EXPECT_FLOAT_EQ(9, stan::math::squared_distance(x2, x1));
  EXPECT_FLOAT_EQ(0, stan::math::squared_distance(x1, x1));
  EXPECT_FLOAT_EQ(0, stan::math::squared_distance(x2, x2));
}

TEST(MathFunctions, squared_distance_nan) {
  double x = 1;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::squared_distance(x, nan), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan, x), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan, nan), std::domain_error);
}

TEST(MathFunctions, squared_distance_inf) {
  double x = 1;
  double inf = std::numeric_limits<double>::infinity();

  EXPECT_THROW(stan::math::squared_distance(x, inf), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf, x), std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf, inf), std::domain_error);
}
