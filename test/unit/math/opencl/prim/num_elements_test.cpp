#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, num_elements_prim) {
  using stan::math::matrix_cl;
  using stan::math::num_elements;
  matrix_cl<double> m0;
  EXPECT_EQ(0, num_elements(m0));

  matrix_cl<double> m1(0, 5);
  EXPECT_EQ(0, num_elements(m1));

  matrix_cl<double> m2(5, 0);
  EXPECT_EQ(0, num_elements(m2));

  matrix_cl<double> m3(5, 4);
  EXPECT_EQ(20, num_elements(m3));

  matrix_cl<double> m4(3, 2);
  EXPECT_EQ(6, num_elements(m4));
}

#endif
