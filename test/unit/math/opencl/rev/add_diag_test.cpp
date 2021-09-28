#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto add_diag_functor
    = [](const auto& a, const auto& b) { return stan::math::add_diag(a, b); };

TEST(OpenCLPrim, add_diag_small) {
  stan::math::matrix_d d1(3, 3);
  d1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::vector_d d2(3);
  d2 << -1, 15, 500;
  stan::math::test::compare_cpu_opencl_prim_rev(add_diag_functor, d1, 2.0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_diag_functor, d1, d2);

  stan::math::matrix_d d3(3, 5);
  d3 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  stan::math::test::compare_cpu_opencl_prim_rev(add_diag_functor, d3, 2);
  stan::math::test::compare_cpu_opencl_prim_rev(add_diag_functor, d3, d2);
}

TEST(OpenCLPrim, add_diag_zero) {
  stan::math::matrix_d d1(0, 0);
  stan::math::vector_d d2(0);
  stan::math::test::compare_cpu_opencl_prim_rev(add_diag_functor, d1, d2);
}

TEST(OpenCLPrim, add_diag_exception) {
  using stan::math::add_diag;
  using stan::math::matrix_cl;
  using stan::math::to_matrix_cl;
  using stan::math::var_value;
  stan::math::matrix_d d1(3, 3);
  stan::math::matrix_d d2(2, 3);
  stan::math::vector_d d3(5);
  stan::math::vector_d d4(3);
  EXPECT_THROW(add_diag(to_matrix_cl(d1), to_matrix_cl(d3)),
               std::invalid_argument);
  EXPECT_THROW(add_diag(to_matrix_cl(d2), to_matrix_cl(d4)),
               std::invalid_argument);
  var_value<matrix_cl<double>> m1_cl = to_matrix_cl(d1);
  var_value<matrix_cl<double>> m2_cl = to_matrix_cl(d2);
  var_value<matrix_cl<double>> m3_cl = to_matrix_cl(d3);
  var_value<matrix_cl<double>> m4_cl = to_matrix_cl(d4);

  EXPECT_THROW(add_diag(m1_cl, m3_cl), std::invalid_argument);
  EXPECT_THROW(add_diag(m2_cl, m4_cl), std::invalid_argument);
}

#endif
