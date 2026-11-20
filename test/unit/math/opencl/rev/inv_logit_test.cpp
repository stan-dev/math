#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto inv_logit_functor = [](const auto& a) { return stan::math::inv_logit(a); };

TEST(OpenCLinv_logit, prim_rev_values_small) {
  Eigen::VectorXd a(14);
  a << -15.2, -10, -0.5, 0.5, 1, 1.0, 1.3, 5, 10, -2.6, -2, -0.2, 1.3, 3;
  stan::math::test::compare_cpu_opencl_prim_rev(inv_logit_functor, a);
}

TEST(OpenCLinv_logit, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(inv_logit_functor, a);
}

TEST(OpenCLinv_logit, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a
      = Eigen::MatrixXd::Constant(N, N, 1.0) + Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(inv_logit_functor, a);
}

TEST(OpenCLinv_logit, rev_derivative_large_args) {
  using stan::math::matrix_cl;
  using stan::math::var_value;
  Eigen::VectorXd x(8);
  x << -40, -20, 0, 20, 30, 37, 40, 100;
  // inv_logit'(x): mpmath, 50 digits
  Eigen::VectorXd d(8);
  d << 4.2483542552915889592e-18, 2.0611536139418493437e-9, 0.25,
      2.0611536139418493437e-9, 9.3576229688384233028e-14,
      8.533047625744064338e-17, 4.2483542552915889592e-18,
      3.720075976020835963e-44;
  const var_value<matrix_cl<double>> x_cl = stan::math::to_matrix_cl(x);
  stan::math::sum(stan::math::inv_logit(x_cl)).grad();
  Eigen::VectorXd adj = stan::math::from_matrix_cl(x_cl.adj());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_NEAR(adj(i), d(i), 1e-13 * d(i)) << "x = " << x(i);
  }
  stan::math::recover_memory();
}

TEST(OpenCLinv_logit, rev_non_const_lvalue) {
  using stan::math::matrix_cl;
  using stan::math::var_value;
  Eigen::VectorXd x(2);
  x << 0, 40;
  var_value<matrix_cl<double>> x_cl = stan::math::to_matrix_cl(x);
  var_value<matrix_cl<double>> y_cl = stan::math::inv_logit(x_cl);
  stan::math::sum(y_cl).grad();
  Eigen::VectorXd adj = stan::math::from_matrix_cl(x_cl.adj());
  EXPECT_NEAR(adj(0), 0.25, 1e-15);
  EXPECT_NEAR(adj(1), 4.2483542552915889592e-18, 1e-30);
  stan::math::recover_memory();
}

#endif
