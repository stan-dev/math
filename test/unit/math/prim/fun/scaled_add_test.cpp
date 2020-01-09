#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, scaled_add) {
  std::vector<double> x(3), y(3);
  double lambda;

  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  y[0] = 2;
  y[1] = 3;
  y[2] = 4;

  lambda = 0.5;

  EXPECT_NO_THROW(stan::math::scaled_add(x, y, lambda));
  EXPECT_FLOAT_EQ(1.0, x[0]);
  EXPECT_FLOAT_EQ(1.5, x[1]);
  EXPECT_FLOAT_EQ(2.0, x[2]);
}

TEST(MathFunctions, scaled_add_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(3), y(3);
  double lambda;

  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  y[0] = 2;
  y[1] = 3;
  y[2] = 4;

  lambda = 0.5;

  EXPECT_NO_THROW(stan::math::scaled_add(x, y, nan));
  EXPECT_TRUE(std::isnan(x[0]));
  EXPECT_TRUE(std::isnan(x[1]));
  EXPECT_TRUE(std::isnan(x[2]));

  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  y[1] = nan;
  EXPECT_NO_THROW(stan::math::scaled_add(x, y, lambda));
  EXPECT_FLOAT_EQ(1.0, x[0]);
  EXPECT_TRUE(std::isnan(x[1]));
  EXPECT_FLOAT_EQ(2.0, x[2]);

  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  EXPECT_NO_THROW(stan::math::scaled_add(x, y, nan));
  EXPECT_TRUE(std::isnan(x[0]));
  EXPECT_TRUE(std::isnan(x[1]));
  EXPECT_TRUE(std::isnan(x[2]));
}
