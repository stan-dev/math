#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto append_col_functor
    = [](const auto& a, const auto& b) { return stan::math::append_col(a, b); };

TEST(OpenCLAppendCol, prim_rev_values_small) {
  int N = 2;
  int M1 = 3;
  int M2 = 2;

  Eigen::MatrixXd a(N, M1);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(N, M2);
  b << 12, 11, 10, 9;
  stan::math::test::compare_cpu_opencl_prim_rev(append_col_functor, a, b);
}

TEST(OpenCLAppendCol, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;

  Eigen::MatrixXd a(N, M);
  Eigen::MatrixXd b(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(append_col_functor, a, b);

  Eigen::MatrixXd c(M, N);
  Eigen::MatrixXd d(M, N);
  stan::math::test::compare_cpu_opencl_prim_rev(append_col_functor, c, d);
}

TEST(OpenCLAppendCol, prim_rev_values_large) {
  int N = 71;
  int M1 = 83;
  int M2 = 93;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M1);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, M2);
  stan::math::test::compare_cpu_opencl_prim_rev(append_col_functor, a, b);
}

#endif
