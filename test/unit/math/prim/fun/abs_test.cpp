#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, abs) {
  using stan::math::abs;

  double y = 2.0;
  EXPECT_FLOAT_EQ(2.0, abs(y));

  y = 128745.72;
  EXPECT_FLOAT_EQ(128745.72, abs(y));

  y = -y;
  EXPECT_FLOAT_EQ(128745.72, abs(y));

  y = -1.3;
  EXPECT_FLOAT_EQ(1.3, abs(y));

  // promoted to double by abs(double)
  int z = 10;
  EXPECT_FLOAT_EQ(10.0, abs(z));
}

TEST(MathFunctions, abs2) {
  double yy = 0;
  yy = 0;
  EXPECT_FLOAT_EQ(0, stan::math::abs(yy));
}

TEST(MathFunctions, abs_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::abs(nan)));
}
