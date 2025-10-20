#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(laplace_marginal_beg_binomial_log_lpmf, phi_dim_2) {
  using stan::math::laplace_marginal_neg_binomial_2_log_lpmf;
  using stan::math::laplace_marginal_tol_neg_binomial_2_log_lpmf;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  double alpha_dbl = 1.6;
  double rho_dbl = 0.45;
  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2);
  x_0 << 0.05100797, 0.16086164;
  Eigen::VectorXd x_1(2);
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  std::vector<double> delta;
  std::vector<int> delta_int;

  std::vector<int> y = {1, 0};
  std::vector<int> y_index = {1, 2};
  double eta_dbl = 10000;

  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& alpha, auto&& rho, auto&& eta) {
          return laplace_marginal_tol_neg_binomial_2_log_lpmf(
              y, y_index, eta, 0, stan::math::test::squared_kernel_functor{},
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(f, alpha_dbl, rho_dbl, eta_dbl);
      },
      theta_0);
}

TEST_F(laplace_disease_map_test, laplace_marginal_neg_binomial_2_log_lpmf) {
  using stan::is_var_v;
  using stan::math::laplace_marginal_neg_binomial_2_log_lpmf;
  using stan::math::laplace_marginal_tol_neg_binomial_2_log_lpmf;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;
  double eta = 1;

  double marginal_density = laplace_marginal_neg_binomial_2_log_lpmf(
      y, y_index, eta, mean, stan::math::test::sqr_exp_kernel_functor(),
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  // ToDo (charlesm93): get benchmark from GPStuff or another software.
  constexpr double tolerance = 1e-6;
  constexpr int max_num_steps = 100;
  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& alpha, auto&& rho, auto&& eta) {
          return laplace_marginal_tol_neg_binomial_2_log_lpmf(
              y, y_index, eta, mean, stan::math::test::sqr_exp_kernel_functor{},
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        auto ret = f(phi_dbl[0], phi_dbl[1], eta);
      },
      theta_0);
  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& alpha, auto&& rho, auto&& eta) {
          return laplace_marginal_tol_neg_binomial_2_log_lpmf(
              y, y_index, eta, mean, stan::math::test::sqr_exp_kernel_functor{},
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(f, phi_dbl[0], phi_dbl[1], eta);
      },
      theta_0);
}
