#ifdef STAN_OPENCL
#include <stan/math/opencl/prim/rows.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::matrix_cl;
using stan::math::rows;
using stan::math::var;

TEST(MathMatrixCL, rows) {
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

  matrix_cl<var> m5(5, 4);
  EXPECT_EQ(5, rows(m5));

  matrix_cl<var> m6(3, 2);
  EXPECT_EQ(3, rows(m6));
}

#endif
