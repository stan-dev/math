#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto append_row_functor
    = [](const auto& a, const auto& b) { return stan::math::append_row(a, b); };

TEST(OpenCLMatrix_append_row, prim_rev_values_small) {
  int N1 = 2;
  int N2 = 3;
  int M = 3;

  Eigen::MatrixXd a(N1, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(N2, M);
  b << 12, 11, 10, 9, 8, 7, 6, 5, 4;
  stan::math::test::compare_cpu_opencl_prim_rev(append_row_functor, a, b);
}

TEST(OpenCLAppendRow, prim_rev_values_M_0) {
  int N1 = 2;
  int N2 = 3;
  int M = 0;

  Eigen::MatrixXd a(N1, M);
  Eigen::MatrixXd b(N2, M);
  stan::math::test::compare_cpu_opencl_prim_rev(append_row_functor, a, b);

  Eigen::MatrixXd c(M, N1);
  Eigen::MatrixXd d(M, N1);
  stan::math::test::compare_cpu_opencl_prim_rev(append_row_functor, c, d);
}

TEST(OpenCLAppendRow, prim_rev_values_large) {
  int N1 = 71;
  int N2 = 95;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N1, M);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N2, M);
  stan::math::test::compare_cpu_opencl_prim_rev(append_row_functor, a, b);
}

#endif
