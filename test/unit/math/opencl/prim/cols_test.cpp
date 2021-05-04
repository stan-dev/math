#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, cols_prim) {
  using stan::math::cols;
  using stan::math::matrix_cl;
  matrix_cl<double> m0;
  EXPECT_EQ(0, cols(m0));

  matrix_cl<double> m1(0, 5);
  EXPECT_EQ(5, cols(m1));

  matrix_cl<double> m2(5, 0);
  EXPECT_EQ(0, cols(m2));

  matrix_cl<double> m3(5, 4);
  EXPECT_EQ(4, cols(m3));

  matrix_cl<double> m4(3, 2);
  EXPECT_EQ(2, cols(m4));
}

#endif
