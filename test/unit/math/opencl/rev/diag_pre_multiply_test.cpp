#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto diag_pre_multiply_functor = [](const auto& a, const auto& b) {
  return stan::math::diag_pre_multiply(a, b);
};

TEST(OpenCL_diag_pre_multiply, diag_pre_multiply_small_vector) {
  Eigen::VectorXd in1(4);
  in1 << 0.5, 3.4, 5.2, 7.5;
  Eigen::MatrixXd in2(4, 2);
  in2 << 3.3, 0.9, 6.7, 1.8, 1, 2, 3, 4;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, in1,
                                                in2);
}

TEST(OpenCL_diag_pre_multiply, diag_pre_multiply_small_row_vector) {
  Eigen::RowVectorXd in1(4);
  in1 << 0.5, 3.4, 5.2, 7.5;
  Eigen::MatrixXd in2(4, 2);
  in2 << 3.3, 0.9, 6.7, 1.8, 1, 2, 3, 4;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, in1,
                                                in2);
}

TEST(OpenCL_diag_pre_multiply, zero) {
  Eigen::VectorXd in1;
  Eigen::MatrixXd in2;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, in1,
                                                in2);
}

TEST(OpenCL_diag_pre_multiply, prim_rev_values_large) {
  int N = 71;

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_pre_multiply_functor, a,
                                                b);
}

#endif
