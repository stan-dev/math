#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, sub_col_exception) {
  stan::math::matrix_cl<double> m1_cl(3, 3);
  EXPECT_THROW(sub_col(m1_cl, 0, 2, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m1_cl, 2, 0, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m1_cl, 4, 2, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m1_cl, 2, 4, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m1_cl, 1, 1, 5), std::invalid_argument);

  stan::math::var_value<stan::math::matrix_cl<double>> m2_cl = m1_cl;
  EXPECT_THROW(sub_col(m2_cl, 0, 2, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m2_cl, 2, 0, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m2_cl, 4, 2, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m2_cl, 2, 4, 1), std::invalid_argument);
  EXPECT_THROW(sub_col(m2_cl, 1, 1, 5), std::invalid_argument);
}

auto sub_col_functor = [](const auto& a, size_t i, size_t j, size_t nrows) {
  return stan::math::sub_col(a, i, j, nrows);
};

TEST(MathMatrixCL, sub_col_value_check) {
  stan::math::matrix_d m1(4, 4);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

  stan::math::test::compare_cpu_opencl_prim_rev(sub_col_functor, m1, 1, 2, 3);
}

#endif
