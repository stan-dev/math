#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, quantile) {
  using stan::math::quantile;

  std::vector<double> v0(3);
  v0[0] = 4;
  v0[1] = 1;
  v0[2] = 2;
  EXPECT_FLOAT_EQ(1.0, quantile(v0, 0.00));
  EXPECT_FLOAT_EQ(1.5, quantile(v0, 0.25));
  EXPECT_FLOAT_EQ(2.0, quantile(v0, 0.50));
  EXPECT_FLOAT_EQ(3.0, quantile(v0, 0.75));
  EXPECT_FLOAT_EQ(4.0, quantile(v0, 1.00));

  std::vector<double> v1(4);
  v1[0] = 7;
  v1[1] = 1;
  v1[2] = 4;
  v1[3] = 2;
  EXPECT_FLOAT_EQ(1.00, quantile(v1, 0.00));
  EXPECT_FLOAT_EQ(1.75, quantile(v1, 0.25));
  EXPECT_FLOAT_EQ(3.00, quantile(v1, 0.50));
  EXPECT_FLOAT_EQ(4.75, quantile(v1, 0.75));
  EXPECT_FLOAT_EQ(7.00, quantile(v1, 1.00));
}
