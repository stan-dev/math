#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto diag_matrix_functor
    = [](const auto& a) { return stan::math::diag_matrix(a).eval(); };

TEST(OpenCLMatrixMultiply, prim_rev_values_small) {
  int N = 2;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_matrix_functor, a);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_N_0) {
  int N = 0;

  Eigen::VectorXd a(N);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_matrix_functor, a);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_large) {
  int N = 71;
  
  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_matrix_functor, a);
}

#endif
