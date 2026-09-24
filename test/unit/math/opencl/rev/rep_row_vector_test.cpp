#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>

auto rep_row_vector_functorCPU
    = [](const auto& a, int n) { return stan::math::rep_row_vector(a, n); };
auto rep_row_vector_functorCL = [](const auto& a, int n) {
  return stan::math::rep_row_vector<stan::conditional_var_value_t<
      decltype(a), stan::math::matrix_cl<double>>>(a, n);
};

TEST(OpenCLRepRowVector, scalar_prim_rev_values_small) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_row_vector_functorCPU, rep_row_vector_functorCL, 6.7, 7);
}

TEST(OpenCLRepRowVector, scalar_prim_rev_size_0) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_row_vector_functorCPU, rep_row_vector_functorCL, 6.7, 0);
}

TEST(OpenCLRepRowVector, scalar_prim_rev_values_large) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_row_vector_functorCPU, rep_row_vector_functorCL, 6.7, 79);
}

#endif
