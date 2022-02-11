#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto to_matrix_functor
    = [](const auto&... a) { return stan::math::to_matrix(a...); };

TEST(OpenCLToMatrix, errors) {
  stan::math::matrix_cl<double> a(2, 3);
  EXPECT_THROW(to_matrix(a, 4, 2), std::invalid_argument);
}

TEST(OpenCLToMatrix, prim_rev_values_one_arg_small) {
  Eigen::VectorXd a(6);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3;
  Eigen::RowVectorXd b = a;
  Eigen::MatrixXd c(3, 2);
  c << 1, 2, 3, 4, 5, 6;
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, b);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, c);
}

TEST(OpenCLToMatrix, prim_rev_values_three_arg_small) {
  Eigen::VectorXd a(6);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3;
  Eigen::RowVectorXd b = a;
  Eigen::MatrixXd c(3, 2);
  c << 1, 2, 3, 4, 5, 6;
  std::vector<double> d{1, -2.3, 3.4, 4, -5.5, 5.7};
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a, 2, 3);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, b, 2, 3);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, c, 2, 3);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, d, 2, 3);
}

TEST(OpenCLToMatrix, prim_rev_values_four_arg_small) {
  std::vector<double> a{1, -2.3, 3.4, 4, -5.5, 5.7};
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a, 2, 3,
                                                false);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a, 2, 3,
                                                true);
}

TEST(OpenCLToMatrix, prim_rev_one_arg_size_0) {
  Eigen::VectorXd a(0);
  Eigen::RowVectorXd b = a;
  Eigen::MatrixXd c(0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, b);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, c);
}

TEST(OpenCLToMatrix, prim_rev_three_arg_size_0) {
  Eigen::VectorXd a(0);
  Eigen::RowVectorXd b = a;
  Eigen::MatrixXd c(0, 0);
  std::vector<double> d{};
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a, 0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, b, 0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, c, 0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, d, 0, 0);
}
TEST(OpenCLToMatrix, prim_rev_four_arg_size_0) {
  std::vector<double> d{};
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, d, 0, 0,
                                                true);
}

TEST(OpenCLToMatrix, prim_rev_values_one_arg_large) {
  int N = 71;
  int M = 87;

  Eigen::VectorXd a(N);
  Eigen::RowVectorXd b(N);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(N, M);

  for (int i = 0; i < N; i++) {
    a[i] = b[i] = Eigen::VectorXd::Random(1)[0];
  }

  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, b);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, c);
}

TEST(OpenCLToMatrix, prim_rev_values_three_arg_large) {
  int N = 71;
  int M = 87;

  Eigen::VectorXd a(N * M);
  Eigen::RowVectorXd b(N * M);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(N, M);
  std::vector<double> d(N * M);

  for (int i = 0; i < N * M; i++) {
    a[i] = b[i] = d[i] = Eigen::VectorXd::Random(1)[0];
  }

  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a, M, N);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, b, M, N);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, c, M, N);
  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, d, M, N);
}

TEST(OpenCLToMatrix, prim_rev_values_four_arg_large) {
  int N = 71;
  int M = 87;

  std::vector<double> a(N * M);

  for (int i = 0; i < N * M; i++) {
    a[i] = Eigen::VectorXd::Random(1)[0];
  }

  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a, M, N,
                                                false);

  stan::math::test::compare_cpu_opencl_prim_rev(to_matrix_functor, a, M, N,
                                                true);
}

#endif
