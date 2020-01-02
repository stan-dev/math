#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, sumZeroSize) {
  std::vector<double> x;
  EXPECT_FLOAT_EQ(0.0, stan::math::sum(x));
}

TEST(MathFunctions, sum) {
  std::vector<double> x(3);

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;

  EXPECT_FLOAT_EQ(6.0, stan::math::sum(x));
}

TEST(MathFunctions, sum_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(3);

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = nan;

  EXPECT_TRUE(std::isnan(stan::math::sum(x)));
}
TEST(MathMatrixPrimArr, sum_vector_int) {
  std::vector<int> x(3);
  EXPECT_EQ(0, stan::math::sum(x));
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  EXPECT_EQ(6, stan::math::sum(x));
}
