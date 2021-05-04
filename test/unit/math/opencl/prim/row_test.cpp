#ifdef STAN_OPENCL
#include <stan/math/opencl/prim/row.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, row_exception) {
  stan::math::matrix_cl<double> m1_cl(3, 3);
  EXPECT_THROW(row(m1_cl, 0), std::invalid_argument);
  EXPECT_THROW(row(m1_cl, 4), std::invalid_argument);
}

TEST(MathMatrixCL, row_value_check) {
  stan::math::matrix_d m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  stan::math::matrix_cl<double> m1_cl(m1);
  stan::math::matrix_cl<double> v_cl = stan::math::row(m1_cl, 3);
  stan::math::matrix_cl<double> v_expr_cl = stan::math::row(m1_cl + m1_cl, 3);

  auto m3 = stan::math::from_matrix_cl(v_cl);
  EXPECT_EQ(7, m3(0));
  EXPECT_EQ(8, m3(1));
  EXPECT_EQ(9, m3(2));

  auto m3_expr = stan::math::from_matrix_cl(v_expr_cl);
  EXPECT_EQ(14, m3_expr(0));
  EXPECT_EQ(16, m3_expr(1));
  EXPECT_EQ(18, m3_expr(2));
}

#endif
