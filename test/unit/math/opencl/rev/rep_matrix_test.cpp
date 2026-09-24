#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>

auto rep_matrix_functorCPU = [](const auto& a, int n, int m) {
  return stan::math::rep_matrix(a, n, m);
};
auto rep_matrix_functorCL = [](const auto& a, int n, int m) {
  return stan::math::rep_matrix<stan::conditional_var_value_t<
      decltype(a), stan::math::matrix_cl<double>>>(a, n, m);
};

auto rep_matrix_functor
    = [](const auto& a, int m) { return stan::math::rep_matrix(a, m); };

TEST(OpenCLRepMatrix, scalar_prim_rev_values_small) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_matrix_functorCPU, rep_matrix_functorCL, 6.7, 7, 3);
}

TEST(OpenCLRepMatrix, scalar_prim_rev_size_0) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_matrix_functorCPU, rep_matrix_functorCL, 6.7, 0, 0);
}

TEST(OpenCLRepMatrix, scalar_prim_rev_values_large) {
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      rep_matrix_functorCPU, rep_matrix_functorCL, 6.7, 79, 83);
}

TEST(OpenCLRepMatrix, vector_prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  Eigen::RowVectorXd b = a;
  stan::math::test::compare_cpu_opencl_prim_rev(rep_matrix_functor, a, 3);
  stan::math::test::compare_cpu_opencl_prim_rev(rep_matrix_functor, b, 3);
}

TEST(OpenCLRepMatrix, vector_prim_rev_size_0) {
  Eigen::VectorXd a(0);
  Eigen::RowVectorXd b = a;
  stan::math::test::compare_cpu_opencl_prim_rev(rep_matrix_functor, a, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(rep_matrix_functor, b, 0);
}

TEST(OpenCLRepMatrix, vector_prim_rev_values_large) {
  int N = 79;
  int M = 85;
  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::RowVectorXd b = a;
  stan::math::test::compare_cpu_opencl_prim_rev(rep_matrix_functor, a, M);
  stan::math::test::compare_cpu_opencl_prim_rev(rep_matrix_functor, b, M);
}

TEST(OpenCLRepMatrix, vector_triangular_adj) {
  using stan::math::matrix_cl;
  using stan::math::var_value;
  Eigen::VectorXd a(8);
  a << -2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4;
  Eigen::MatrixXd adj = Eigen::MatrixXd::Ones(8, 8);
  Eigen::RowVectorXd b = a;
  matrix_cl<double> a_cl(a);
  matrix_cl<double> b_cl(b);
  matrix_cl<double> adj_cl(adj);
  var_value<matrix_cl<double>> a_var(a_cl);
  var_value<matrix_cl<double>> b_var(b_cl);

  var_value<matrix_cl<double>> a_res = stan::math::rep_matrix(a_var, 8);
  var_value<matrix_cl<double>> b_res = stan::math::rep_matrix(b_var, 8);
  adj_cl.view(stan::math::matrix_cl_view::Lower);
  a_res.adj() = adj_cl;
  b_res.adj() = adj_cl;
  stan::math::grad();
  EXPECT_MATRIX_EQ(from_matrix_cl(a_var.adj()),
                   stan::math::linspaced_vector(8, 1., 8.));
  EXPECT_MATRIX_EQ(from_matrix_cl(b_var.adj()),
                   stan::math::linspaced_row_vector(8, 1., 8.).reverse());

  stan::math::set_zero_all_adjoints();
  adj_cl.view(stan::math::matrix_cl_view::Upper);
  a_res.adj() = adj_cl;
  b_res.adj() = adj_cl;
  stan::math::grad();
  EXPECT_MATRIX_EQ(from_matrix_cl(a_var.adj()),
                   stan::math::linspaced_vector(8, 1., 8.).reverse());
  EXPECT_MATRIX_EQ(from_matrix_cl(b_var.adj()),
                   stan::math::linspaced_row_vector(8, 1., 8.));
}

#endif
