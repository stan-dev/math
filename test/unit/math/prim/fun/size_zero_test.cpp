#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixPrimScal, size_zero) {
  std::vector<double> x(3), y(0);
  stan::math::matrix_d m(3, 4);
  stan::math::matrix_d n(3, 0);
  int f, g;

  EXPECT_EQ(0, stan::math::size_zero(x));
  EXPECT_EQ(0, stan::math::size_zero(m));
  EXPECT_EQ(0, stan::math::size_zero(x, f, g, m));

  EXPECT_EQ(1, stan::math::size_zero(y));
  EXPECT_EQ(1, stan::math::size_zero(n));
  EXPECT_EQ(1, stan::math::size_zero(y, n));
  EXPECT_EQ(1, stan::math::size_zero(x, f, y, g));
}
