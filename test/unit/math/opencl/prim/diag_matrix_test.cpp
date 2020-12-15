#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto diag_matrix_functor
    = [](const auto& a) { return stan::math::diag_matrix(a); };

TEST(OpenCLPrim, diag_matrix_small) {
  stan::math::vector_d d3(3);
  d3 << 1, 2, 3;
  stan::math::test::compare_cpu_opencl_prim(diag_matrix_functor, d3);

  stan::math::vector_d d0(0);
  stan::math::test::compare_cpu_opencl_prim(diag_matrix_functor, d0);

  stan::math::vector_d d1(1);
  d1 << 12;
  stan::math::test::compare_cpu_opencl_prim(diag_matrix_functor, d1);
}

#endif
