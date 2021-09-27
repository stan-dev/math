#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto diag_post_multiply_functor = [](const auto& a, const auto& b) {
  return stan::math::diag_post_multiply(a, b);
};

TEST(OpenCL_diag_post_multiply, diag_post_multiply_small_vector) {
  Eigen::MatrixXd in1(4, 2);
  in1 << 3.3, 0.9, 6.7, 1.8, 1, 2, 3, 4;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, in1,
                                                in2);
}

TEST(OpenCL_diag_post_multiply, diag_post_multiply_small_row_vector) {
  Eigen::MatrixXd in1(4, 2);
  in1 << 3.3, 0.9, 6.7, 1.8, 1, 2, 3, 4;
  Eigen::RowVectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, in1,
                                                in2);
}

TEST(OpenCL_diag_post_multiply, zero) {
  Eigen::MatrixXd in1(4, 0);
  Eigen::VectorXd in2;
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, in1,
                                                in2);
}

TEST(OpenCL_diag_post_multiply, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  Eigen::RowVectorXd b = Eigen::RowVectorXd::Random(N);
  stan::math::test::compare_cpu_opencl_prim_rev(diag_post_multiply_functor, a,
                                                b);
}

#endif
