#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, rows_prim) {
  using stan::math::matrix_cl;
  using stan::math::rows;
  matrix_cl<double> m0;
  EXPECT_EQ(0, rows(m0));

  matrix_cl<double> m1(0, 5);
  EXPECT_EQ(0, rows(m1));

  matrix_cl<double> m2(5, 0);
  EXPECT_EQ(5, rows(m2));

  matrix_cl<double> m3(5, 4);
  EXPECT_EQ(5, rows(m3));

  matrix_cl<double> m4(3, 2);
  EXPECT_EQ(3, rows(m4));
}

#endif
