#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, quantile) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::quantile;

  Matrix<double, Dynamic, Dynamic> m(3, 2);
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_FLOAT_EQ(1.0, quantile(m, 0.0));
  EXPECT_FLOAT_EQ(3.5, quantile(m, 0.5));
  EXPECT_FLOAT_EQ(6.0, quantile(m, 1.0));

  Matrix<double, Dynamic, 1> v0(3);
  v0 << 4, 1, 2;
  EXPECT_FLOAT_EQ(1.0, quantile(v0, 0.00));
  EXPECT_FLOAT_EQ(1.5, quantile(v0, 0.25));
  EXPECT_FLOAT_EQ(2.0, quantile(v0, 0.50));
  EXPECT_FLOAT_EQ(3.0, quantile(v0, 0.75));
  EXPECT_FLOAT_EQ(4.0, quantile(v0, 1.00));

  Matrix<double, Dynamic, 1> v1(4);
  v1 << 7, 1, 4, 2;
  EXPECT_FLOAT_EQ(1.00, quantile(v1, 0.00));
  EXPECT_FLOAT_EQ(1.75, quantile(v1, 0.25));
  EXPECT_FLOAT_EQ(3.00, quantile(v1, 0.50));
  EXPECT_FLOAT_EQ(4.75, quantile(v1, 0.75));
  EXPECT_FLOAT_EQ(7.00, quantile(v1, 1.00));
}
