#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto to_row_vector_functor
    = [](const auto& a) { return stan::math::to_row_vector(a); };

TEST(OpenCLToRowVector, prim_rev_values_small) {
  Eigen::VectorXd a(6);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3;
  Eigen::RowVectorXd b = a;
  Eigen::MatrixXd c(3, 2);
  c << 1, 2, 3, 4, 5, 6;
  std::vector<double> d{1, -2.3, 3.4, 4, -5.5, 5.7};
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, b);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, c);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, d);
}

TEST(OpenCLToRowVector, prim_rev_size_0) {
  Eigen::VectorXd a(0);
  Eigen::RowVectorXd b = a;
  Eigen::MatrixXd c(0, 0);
  std::vector<double> d{};
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, b);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, c);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, d);
}

TEST(OpenCLToRowVector, prim_rev_values_large) {
  int N = 71;
  int M = 87;

  Eigen::VectorXd a(N);
  Eigen::RowVectorXd b(N);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(N, M);
  std::vector<double> d(N * M);

  for (int i = 0; i < N; i++) {
    a[i] = b[i] = d[i] = Eigen::VectorXd::Random(1)[0];
  }

  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, a);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, b);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, c);
  stan::math::test::compare_cpu_opencl_prim_rev(to_row_vector_functor, d);
}

#endif
