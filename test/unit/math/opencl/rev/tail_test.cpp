#ifdef STAN_OPENCL
#include <stan/math/rev.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixTailCL, tail_size) {
  using stan::math::matrix_cl;
  using stan::math::tail;
  using stan::math::var_value;
  matrix_cl<double> v(3, 1);
  EXPECT_EQ(0, tail(v, 0).eval().size());
  EXPECT_EQ(3, tail(v, 3).eval().size());
  EXPECT_THROW(tail(v, 4), std::out_of_range);

  matrix_cl<double> rv(1, 3);
  EXPECT_EQ(0, tail(rv, 0).eval().size());
  EXPECT_EQ(3, tail(rv, 3).eval().size());
  EXPECT_THROW(tail(rv, 4), std::out_of_range);

  var_value<matrix_cl<double>> v_var = v;
  EXPECT_EQ(0, tail(v_var, 0).size());
  EXPECT_EQ(3, tail(v_var, 3).size());
  EXPECT_THROW(tail(v_var, 4), std::out_of_range);

  var_value<matrix_cl<double>> rv_var = rv;
  EXPECT_EQ(0, tail(rv_var, 0).size());
  EXPECT_EQ(3, tail(rv_var, 3).size());
  EXPECT_THROW(tail(rv_var, 4), std::out_of_range);
}

auto tail_functor = [](const auto& a) { return stan::math::tail(a, 5); };

TEST(MathMatrixTailCL, tail_value_check_vector) {
  stan::math::vector_d m1(9);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  stan::math::test::compare_cpu_opencl_prim_rev(tail_functor, m1);
}

TEST(MathMatrixTailCL, tail_value_check_row_vector) {
  stan::math::row_vector_d m1(9);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  stan::math::test::compare_cpu_opencl_prim_rev(tail_functor, m1);
}
#endif
