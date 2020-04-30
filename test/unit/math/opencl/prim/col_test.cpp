#ifdef STAN_OPENCL
#include <stan/math/opencl/prim/col.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, col_exception) {
  stan::math::matrix_cl<double> m1_cl(3,3);
  EXPECT_THROW(col(m1_cl, 0), std::invalid_argument);
  EXPECT_THROW(col(m1_cl, 4), std::domain_error);
}


TEST(MathMatrixCL, col_value_check) {
  stan::math::matrix_d m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  stan::math::matrix_cl<double> m1_cl(m1);
  stan::math::matrix_cl<double> v_cl = stan::math::col(m1_cl, 3);
  
  stan::math::vector_d m3 = stan::math::from_matrix_cl(v_cl);
  EXPECT_EQ(3, m3(0, 0));
  EXPECT_EQ(6, m3(1, 0));
  EXPECT_EQ(9, m3(2, 0));
}

#endif
