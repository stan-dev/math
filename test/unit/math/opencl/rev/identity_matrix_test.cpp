#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>

auto identity_matrix_functorCPU
    = [](int n) { return stan::math::identity_matrix(n); };
auto identity_matrix_functorCL = [](int n) {
  return stan::math::identity_matrix<stan::math::matrix_cl<double>>(n);
};

TEST(OpenCLIdentityMatrix, scalar_prim_rev_values_small) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      identity_matrix_functorCPU, identity_matrix_functorCL, 7);
}

TEST(OpenCLIdentityMatrix, scalar_prim_rev_size_0) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      identity_matrix_functorCPU, identity_matrix_functorCL, 0);
}

TEST(OpenCLIdentityMatrix, scalar_prim_rev_values_large) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      identity_matrix_functorCPU, identity_matrix_functorCL, 79);
}

#endif
