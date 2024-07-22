#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix/laplace.hpp>
#include <test/unit/math/mix/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace_marginal_poisson_log_lpmf, phi_dim_2) {
  using stan::math::laplace_marginal_poisson_2_log_lpmf;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_2_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;

  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  int dim_phi = 2;
  double alpha_dbl = 1.6;
  double rho_dbl = 0.45;
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

  stan::math::test::squared_kernel_functor sq_kernel;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 4e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  double tolerance = 1e-6;
  int max_num_steps = 100;
  for (int max_steps_line_search = 0; max_steps_line_search < 4; ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3; hessian_block_size++) {
        for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
            return laplace_marginal_tol_poisson_log_lpmf(
            sums, n_samples, tolerance, max_num_steps, hessian_block_size,
            solver_num, max_steps_line_search, theta_0, sq_kernel, nullptr, x, alpha, rho);
        };
        stan::test::expect_ad<true>(ad_tol, f, alpha_dbl, rho_dbl);
        }
    }
  }
  Eigen::VectorXd ye(2);
  ye << 1, 1;
  for (int max_steps_line_search = 0; max_steps_line_search < 4; ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3; hessian_block_size++) {
        for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
            return laplace_marginal_tol_poisson_2_log_lpmf(
            sums, n_samples, ye, tolerance, max_num_steps, hessian_block_size,
            solver_num, max_steps_line_search, theta_0, sq_kernel, nullptr, x, alpha, rho);
        };
        stan::test::expect_ad<true>(ad_tol, f, alpha_dbl, rho_dbl);
        }
    }
  }
}

TEST_F(laplace_disease_map_test, laplace_marginal_poisson_log_lpmf) {
  using stan::math::laplace_marginal_poisson_2_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_2_log_lpmf;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::value_of;
  using stan::math::var;

  double marginal_density = laplace_marginal_poisson_2_log_lpmf(
      y, n_samples, ye, theta_0, stan::math::test::sqr_exp_kernel_functor(),
      nullptr, x, phi_dbl(0), phi_dbl(1));

  double tol = 6e-4;
  // Benchmark from GPStuff.
  EXPECT_NEAR(-2866.88, marginal_density, tol);

  stan::math::test::squared_kernel_functor sq_kernel;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 4e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  double tolerance = 1e-6;
  int max_num_steps = 100;
  for (int max_steps_line_search = 0; max_steps_line_search < 4; ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3; hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
            return laplace_marginal_tol_poisson_2_log_lpmf(
            y, n_samples, ye, tolerance, max_num_steps, hessian_block_size,
            solver_num, max_steps_line_search, theta_0,
            stan::math::test::sqr_exp_kernel_functor(), nullptr, x, alpha, rho);
        };
        stan::test::expect_ad<true>(ad_tol, f, phi_dbl[0], phi_dbl[1]);
      }
    }
  }
}
/*
TEST(laplace_marginal_poisson_log_lpmf, phi_dim_2) {
  using stan::math::laplace_marginal_poisson_2_log_lpmf;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_2_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;

  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  int dim_phi = 2;
  var alpha = 1.6;
  var rho = 0.45;
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

  stan::math::test::squared_kernel_functor sq_kernel;
  var target = laplace_marginal_poisson_log_lpmf(
      sums, n_samples, theta_0, sq_kernel, nullptr, x, alpha, rho);

  // TODO: benchmark target against gpstuff.
  // Expected: -2.53056
  double tol = 1e-4;
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with exposure argument
  Eigen::VectorXd ye(2);
  ye << 1, 1;
  target = laplace_marginal_poisson_2_log_lpmf(
      sums, n_samples, ye, theta_0, sq_kernel, nullptr, x, alpha, rho);
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with optional arguments
  double tolerance = 1e-6;
  int max_num_steps = 100;
  target = laplace_marginal_tol_poisson_2_log_lpmf(
      sums, n_samples, ye, tolerance, max_num_steps, 1, 1, 0, theta_0,
      sq_kernel, nullptr, x, alpha, rho);
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  std::vector<double> g;
  std::vector<stan::math::var> parm_vec{alpha, rho};
  target.grad(parm_vec, g);

  // finite diff test
  double diff = 1e-7;
  Eigen::VectorXd phi_dbl(2);  // = value_of(phi);
  phi_dbl(0) = value_of(alpha);
  phi_dbl(1) = value_of(rho);
  Eigen::VectorXd phi_1l = phi_dbl;
  Eigen::VectorXd phi_1u = phi_dbl;
  Eigen::VectorXd phi_2l = phi_dbl;
  Eigen::VectorXd phi_2u = phi_dbl;
  phi_1l(0) -= diff;
  phi_1u(0) += diff;
  phi_2l(1) -= diff;
  phi_2u(1) += diff;

  double target_1u = laplace_marginal_tol_poisson_log_lpmf(
      sums, n_samples, tolerance, max_num_steps, 1, 1, 0, theta_0, sq_kernel,
      nullptr, x, phi_1u[0], phi_1u[1]);
  double target_1l = laplace_marginal_tol_poisson_log_lpmf(
      sums, n_samples, tolerance, max_num_steps, 1, 1, 0, theta_0, sq_kernel,
      nullptr, x, phi_1l[0], phi_1l[1]);
  double target_2u = laplace_marginal_tol_poisson_log_lpmf(
      sums, n_samples, tolerance, max_num_steps, 1, 1, 0, theta_0, sq_kernel,
      nullptr, x, phi_2u[0], phi_2u[1]);
  double target_2l = laplace_marginal_tol_poisson_log_lpmf(
      sums, n_samples, tolerance, max_num_steps, 1, 1, 0, theta_0, sq_kernel,
      nullptr, x, phi_2l[0], phi_2l[1]);

  std::vector<double> g_finite(dim_phi);
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);

  tol = 1.1e-4;
  EXPECT_NEAR(g_finite[0], g[0], tol);
  EXPECT_NEAR(g_finite[1], g[1], tol);
}
*/
