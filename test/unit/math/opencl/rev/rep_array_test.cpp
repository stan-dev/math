#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>

auto rep_array_functorCPU
    = [](const auto& a, int n) { return stan::math::rep_array(a, n); };
auto rep_array_functorCL = [](const auto& a, int n) {
  return stan::math::rep_array<stan::conditional_var_value_t<
      decltype(a), stan::math::matrix_cl<double>>>(a, n);
};

TEST(OpenCLRepArray, scalar_prim_rev_values_small) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_array_functorCPU, rep_array_functorCL, 6.7, 7);
}

TEST(OpenCLRepArray, scalar_prim_rev_size_0) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_array_functorCPU, rep_array_functorCL, 6.7, 0);
}

TEST(OpenCLRepArray, scalar_prim_rev_values_large) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_array_functorCPU, rep_array_functorCL, 6.7, 79);
}

#endif
