#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto log_inv_logit_diff_functor = [](const auto& a, const auto& b) {
  return stan::math::log_inv_logit_diff(a, b);
};

TEST(OpenCLMatrix_log_inv_logit_diff, prim_rev_values_small) {
  int N = 2;
  int M = 3;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(N, M);
  b << 12, 11, 10, 9, 8, 7;
  stan::math::test::compare_cpu_opencl_prim_rev(log_inv_logit_diff_functor, a,
                                                b);
}

TEST(OpenCLMatrix_log_inv_logit_diff, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;

  Eigen::MatrixXd a(N, M);
  Eigen::MatrixXd b(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(log_inv_logit_diff_functor, a,
                                                b);

  Eigen::MatrixXd c(M, N);
  Eigen::MatrixXd d(M, N);
  stan::math::test::compare_cpu_opencl_prim_rev(log_inv_logit_diff_functor, c,
                                                d);
}

TEST(OpenCLMatrix_log_inv_logit_diff, prim_rev_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(N, M);
  stan::math::test::compare_cpu_opencl_prim_rev(log_inv_logit_diff_functor, a,
                                                b);
}

TEST(OpenCLMatrix_log_inv_logit_diff, prim_rev_scalar_values_large) {
  int N = 71;
  int M = 83;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  double b = 0.3;
  stan::math::test::compare_cpu_opencl_prim_rev(log_inv_logit_diff_functor, a,
                                                b);
  stan::math::test::compare_cpu_opencl_prim_rev(log_inv_logit_diff_functor, b,
                                                a);
}

TEST(OpenCLMatrix_log_inv_logit_diff, rev_gradient_large_args) {
  using stan::math::matrix_cl;
  using stan::math::var_value;
  Eigen::VectorXd x(8);
  Eigen::VectorXd y(8);
  x << 40, 100, 30, 37, -30, 1, 700, 0.5;
  y << 0, 1, -30, 36, -40, 0, 20, -1;
  // d/dx and d/dy: mpmath, 50 digits
  Eigen::VectorXd dx(8);
  Eigen::VectorXd dy(8);
  dx << 8.4967085105831779907e-18, 1.3832290902125321262e-43,
      9.3576229688401746049e-14, 5.8197670686932650972e-1,
      1.0000454019909161115, 8.5091812823932154513e-1,
      4.7835719068902113962e-296, 6.647575855870136797e-1;
  dy << -5.0000000000000000425e-1, -7.3105857863000487925e-1,
      -9.3576229688401746049e-14, -1.5819767068693261924,
      -4.5401991009692016683e-5, -1.0819767068693264244,
      -9.9999999793884638181e-1, -5.5615833815886336509e-1;
  var_value<matrix_cl<double>> x_cl = stan::math::to_matrix_cl(x);
  var_value<matrix_cl<double>> y_cl = stan::math::to_matrix_cl(y);
  stan::math::sum(stan::math::log_inv_logit_diff(x_cl, y_cl)).grad();
  Eigen::VectorXd x_adj = stan::math::from_matrix_cl(x_cl.adj());
  Eigen::VectorXd y_adj = stan::math::from_matrix_cl(y_cl.adj());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_NEAR(x_adj(i), dx(i), 1e-13 * dx(i)) << x(i) << ", " << y(i);
    EXPECT_NEAR(y_adj(i), dy(i), -1e-13 * dy(i)) << x(i) << ", " << y(i);
  }
  stan::math::recover_memory();
}

#endif
