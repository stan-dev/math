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

namespace stan {
namespace test {
template <typename T1, typename T2, typename T3, typename T4, typename T5>
inline void test_laplace_marginal_poisson_logit_lpmf_finite_diff(
    T1&& y, T2&& phi_dbl, T3&& x, T4&& n_samples, T5&& theta_0,
    double tolerance, long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search, const double gp_tol,
    const double fd_tol) {
  using stan::math::diff_bernoulli_logit;
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;
  using stan::math::laplace_marginal_tol_bernoulli_logit_lpmf;
  using stan::math::to_vector;
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi = phi_dbl;
  // Test with optional arguments.
  stan::math::test::sqr_exp_kernel_functor kernel_functor;
  var target = laplace_marginal_tol_bernoulli_logit_lpmf(
      y, n_samples, tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, theta_0, kernel_functor, 0, x, phi(0), phi(1));
  std::string err_msg
      = "laplace parameters: \n"
        "\ttolerance("
        + std::to_string(tolerance) + ")" + "\n\tmax_num_steps("
        + std::to_string(max_num_steps) + ")" + "\n\thessian_block_size("
        + std::to_string(hessian_block_size) + ")" + "\n\tsolver("
        + std::to_string(solver) + ")" + "\n\tmax_steps_line_search("
        + std::to_string(max_steps_line_search) + ")";
  EXPECT_NEAR(-195.368, value_of(target), gp_tol) << err_msg;

  std::vector<double> g;
  std::vector<stan::math::var> parm_vec{phi(0), phi(1)};
  target.grad(parm_vec, g);

  // finite diff benchmark
  double diff = 1e-7;
  Eigen::VectorXd phi_1l = phi_dbl;
  phi_1l(0) -= diff;
  Eigen::VectorXd phi_1u = phi_dbl;
  phi_1u(0) += diff;
  Eigen::VectorXd phi_2l = phi_dbl;
  phi_2l(1) -= diff;
  Eigen::VectorXd phi_2u = phi_dbl;
  phi_2u(1) += diff;

  double target_1u = laplace_marginal_tol_bernoulli_logit_lpmf(
      y, n_samples, tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, theta_0, kernel_functor, nullptr, x, phi_1u(0),
      phi_1u(1));
  double target_1l = laplace_marginal_tol_bernoulli_logit_lpmf(
      y, n_samples, tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, theta_0, kernel_functor, nullptr, x, phi_1l(0),
      phi_1l(1));
  double target_2u = laplace_marginal_tol_bernoulli_logit_lpmf(
      y, n_samples, tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, theta_0, kernel_functor, nullptr, x, phi_2u(0),
      phi_2u(1));
  double target_2l = laplace_marginal_tol_bernoulli_logit_lpmf(
      y, n_samples, tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, theta_0, kernel_functor, nullptr, x, phi_2l(0),
      phi_2l(1));

  std::vector<double> g_finite(phi_dbl.size());
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);
  EXPECT_NEAR(g_finite[0], g[0], fd_tol) << err_msg;
  EXPECT_NEAR(g_finite[1], g[1], fd_tol) << err_msg;
  stan::math::recover_memory();
}
}  // namespace test
}  // namespace stan
TEST(laplace_marginal_bernoulli_logit_lpmf, phi_dim500) {
  using stan::math::diff_bernoulli_logit;
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;
  using stan::math::laplace_marginal_tol_bernoulli_logit_lpmf;
  using stan::math::to_vector;
  using stan::math::var;

  int dim_theta = 500;
  int n_observations = 500;
  std::string data_directory = "test/unit/math/mix/laplace/aki_synth_data/";
  std::vector<double> x1(dim_theta), x2(dim_theta);
  std::vector<int> y(n_observations);
  stan::math::test::read_in_data(dim_theta, n_observations, data_directory, x1,
                                 x2, y);

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    Eigen::VectorXd coordinate(dim_x);
    coordinate << x1[i], x2[i];
    x[i] = coordinate;
  }
  std::vector<int> n_samples = stan::math::rep_array(1, dim_theta);
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  std::vector<double> delta;
  std::vector<int> delta_int;
  int dim_phi = 2;
  double tol = 8e-5;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl(dim_phi);
  phi_dbl << 1.6, 1;
  stan::math::test::sqr_exp_kernel_functor kernel_functor;
  double target = laplace_marginal_bernoulli_logit_lpmf(
      y, n_samples, theta_0, kernel_functor, nullptr, x, phi_dbl(0),
      phi_dbl(1));
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, target, tol);
  double tolerance = 1e-6;
  int max_num_steps = 100;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 4e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  // FIXME(Steve): hessian_block_size of 3 fails approx test
  for (int max_steps_line_search = 0; max_steps_line_search < 4;
       ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3;
         hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_bernoulli_logit_lpmf(
              y, n_samples, tolerance, max_num_steps, hessian_block_size,
              solver_num, max_steps_line_search, theta_0, kernel_functor, 0, x,
              alpha, rho);
        };
        stan::test::expect_ad<true>(ad_tol, f, phi_dbl[0], phi_dbl[1]);
      }
    }
  }
}
