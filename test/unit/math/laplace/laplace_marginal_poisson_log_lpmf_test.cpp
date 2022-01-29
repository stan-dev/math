#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace_marginal_poisson_log_lpmf, phi_dim_2) {
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  int dim_phi = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
  phi << 1.6, 0.45;

  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2);
  x_0 << 0.05100797, 0.16086164;
  Eigen::VectorXd x_1(2);
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  std::vector<double> delta;
  std::vector<int> delta_int;

  std::vector<int> sums = {1, 0};
  std::vector<int> n_samples = {1, 1};

  stan::math::test::squared_kernel_functor K;
  var target = laplace_marginal_poisson_log_lpmf(sums, n_samples, K, phi, x,
                                                 delta, delta_int, theta_0);

  // TODO: benchmark target against gpstuff.
  // Expected: -2.53056
  double tol = 1e-4;
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with exposure argument
  Eigen::VectorXd ye(2);
  ye << 1, 1;
  target = laplace_marginal_poisson_log_lpmf(sums, n_samples, ye, K, phi, x,
                                             delta, delta_int, theta_0);
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with optional arguments
  double tolerance = 1e-6;
  int max_num_steps = 100;
  target = laplace_marginal_poisson_log_lpmf(sums, n_samples, ye, K, phi, x,
                                             delta, delta_int, theta_0,
                                             0, tolerance, max_num_steps);
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  std::vector<double> g;
  std::vector<stan::math::var> parm_vec{phi(0), phi(1)};
  target.grad(parm_vec, g);

  // finite diff test
  double diff = 1e-7;
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_1l = phi_dbl, phi_1u = phi_dbl, phi_2l = phi_dbl,
                  phi_2u = phi_dbl;
  phi_1l(0) -= diff;
  phi_1u(0) += diff;
  phi_2l(1) -= diff;
  phi_2u(1) += diff;

  double target_1u = laplace_marginal_poisson_log_lpmf(
             sums, n_samples, K, phi_1u, x, delta, delta_int, theta_0),
         target_1l = laplace_marginal_poisson_log_lpmf(
             sums, n_samples, K, phi_1l, x, delta, delta_int, theta_0),
         target_2u = laplace_marginal_poisson_log_lpmf(
             sums, n_samples, K, phi_2u, x, delta, delta_int, theta_0),
         target_2l = laplace_marginal_poisson_log_lpmf(
             sums, n_samples, K, phi_2l, x, delta, delta_int, theta_0);

  std::vector<double> g_finite(dim_phi);
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);

  tol = 1.1e-4;
  EXPECT_NEAR(g_finite[0], g[0], tol);
  EXPECT_NEAR(g_finite[1], g[1], tol);
}


TEST_F(laplace_disease_map_test, laplace_marginal_poisson_log_lpmf) {
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::var;
  using stan::math::value_of;

  var marginal_density = laplace_marginal_poisson_log_lpmf(
      y, n_samples, ye, stan::math::test::sqr_exp_kernel_functor(), phi, x,
      delta, delta_int, theta_0);

  double tol = 6e-4;
  // Benchmark from GPStuff.
  EXPECT_NEAR(-2866.88, value_of(marginal_density), tol);

  std::vector<double> g;
  std::vector<var> parm_vec{phi(0), phi(1)};
  marginal_density.grad(parm_vec, g);

  // finite diff
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_u0 = phi_dbl, phi_u1 = phi_dbl, phi_l0 = phi_dbl,
                  phi_l1 = phi_dbl;
  double eps = 1e-7;

  phi_u0(0) += eps;
  phi_u1(1) += eps;
  phi_l0(0) -= eps;
  phi_l1(1) -= eps;

  double target_u0 = laplace_marginal_poisson_log_lpmf(
             y, n_samples, ye, stan::math::test::sqr_exp_kernel_functor(),
             phi_u0, x, delta, delta_int, theta_0),

         target_u1 = laplace_marginal_poisson_log_lpmf(
             y, n_samples, ye, stan::math::test::sqr_exp_kernel_functor(),
             phi_u1, x, delta, delta_int, theta_0),

         target_l0 = laplace_marginal_poisson_log_lpmf(
             y, n_samples, ye, stan::math::test::sqr_exp_kernel_functor(),
             phi_l0, x, delta, delta_int, theta_0),

         target_l1 = laplace_marginal_poisson_log_lpmf(
             y, n_samples, ye, stan::math::test::sqr_exp_kernel_functor(),
             phi_l1, x, delta, delta_int, theta_0);

  EXPECT_NEAR((target_u0 - target_l0) / (2.0 * eps), g[0], 8e-4);
  EXPECT_NEAR((target_u1 - target_l1) / (2.0 * eps), g[1], 0.00017);
}
