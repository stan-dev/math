#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
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
  for (int max_steps_line_search = 0; max_steps_line_search < 4;
       ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3;
         hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_poisson_log_lpmf(
              sums, n_samples, theta_0, sq_kernel,
              std::forward_as_tuple(x, alpha, rho), tolerance, max_num_steps,
              hessian_block_size, solver_num, max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(ad_tol, f, alpha_dbl, rho_dbl);
      }
    }
  }

  Eigen::VectorXd ye(2);
  ye << 1, 1;
  for (int max_steps_line_search = 0; max_steps_line_search < 4;
       ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3;
         hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_poisson_2_log_lpmf(
              sums, n_samples, ye, theta_0, sq_kernel,
              std::forward_as_tuple(x, alpha, rho), tolerance, max_num_steps,
              hessian_block_size, solver_num, max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(ad_tol, f, alpha_dbl, rho_dbl);
      }
    }
  }
}

TEST_F(laplace_disease_map_test, laplace_marginal_poisson_log_lpmf) {
  using stan::math::laplace_marginal_poisson_2_log_lpmf;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_2_log_lpmf;
  using stan::math::value_of;
  using stan::math::var;

  double marginal_density = laplace_marginal_poisson_2_log_lpmf(
      y, n_samples, ye, theta_0, stan::math::test::sqr_exp_kernel_functor(),
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  double tol = 6e-4;
  // Benchmark from GPStuff.
  EXPECT_NEAR(-2866.88, marginal_density, tol);

  stan::math::test::squared_kernel_functor sq_kernel;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 4e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  double tolerance = 1e-6;
  int max_num_steps = 100;
  for (int max_steps_line_search = 0; max_steps_line_search < 4;
       ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3;
         hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_poisson_2_log_lpmf(
              y, n_samples, ye, theta_0,
              stan::math::test::sqr_exp_kernel_functor(),
              std::forward_as_tuple(x, alpha, rho), tolerance, max_num_steps,
              hessian_block_size, solver_num, max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(ad_tol, f, phi_dbl[0], phi_dbl[1]);
      }
    }
  }
}
