#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto rows_dot_product_functor = [](const auto& a, const auto& b) {
  return stan::math::rows_dot_product(a, b);
};

TEST(OpenCLRowsDotProduct, errors) {
  Eigen::MatrixXd a(3, 3);
  Eigen::MatrixXd b(8, 3);
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  EXPECT_THROW(stan::math::rows_dot_product(a_cl, b_cl), std::invalid_argument);
  EXPECT_THROW(stan::math::rows_dot_product(a_cl, b_cl), std::invalid_argument);
}

TEST(OpenCLRowsDotProduct, prim_rev_small_vector) {
  Eigen::MatrixXd a(2, 4);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  Eigen::MatrixXd b(2, 4);
  b << 1, 2, 3, 4, 5, 6, 7, 8;
  stan::math::test::compare_cpu_opencl_prim_rev(rows_dot_product_functor, a, b);
}

TEST(OpenCLRowsDotProduct, prim_rev_size_0) {
  Eigen::MatrixXd a(0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(rows_dot_product_functor, a, a);
}

TEST(OpenCLRowsDotProduct, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(rows_dot_product_functor, a, b);
}

#endif
