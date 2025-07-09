#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/laplace/aki_synth_data/x1.hpp>

#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(laplace_marginal_bernoulli_logit_lpmf, phi_dim500) {
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;
  using stan::math::laplace_marginal_tol_bernoulli_logit_lpmf;
  using stan::math::to_vector;
  using stan::math::var;
  using stan::math::test::flag_test;
  int dim_theta = 500;
  int n_observations = 500;
  auto x1 = stan::test::laplace::x1;
  auto x2 = stan::test::laplace::x2;
  auto y = stan::test::laplace::y;

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    Eigen::VectorXd coordinate(dim_x);
    coordinate << x1[i], x2[i];
    x[i] = coordinate;
  }
  std::vector<int> n_samples = stan::math::rep_array(1, dim_theta);
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  Eigen::VectorXd mean = Eigen::VectorXd::Zero(dim_theta);
  std::vector<double> delta;
  std::vector<int> delta_int;
  int dim_phi = 2;
  double tol = 8e-5;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl(dim_phi);
  phi_dbl << 1.6, 1;
  using stan::math::test::sqr_exp_kernel_functor;
  double target = laplace_marginal_bernoulli_logit_lpmf(
      y, n_samples, 0, sqr_exp_kernel_functor{},
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, target, tol);
  constexpr double tolerance = 1e-8;
  constexpr int max_num_steps = 1000;
  // solver_num, max_steps_line_search, hessian_block_size
  using stan::math::test::laplace_issue;
  constexpr std::array known_issues{
      laplace_issue{1, 1, 1}, laplace_issue{1, 2, 1}, laplace_issue{1, 2, 3}};
  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_bernoulli_logit_lpmf(
              y, n_samples, mean, sqr_exp_kernel_functor{},
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::ad_tolerances tols;
        if (flag_test(known_issues, solver_num, max_steps_line_search,
                      hessian_block_size)) {
          tols.gradient_grad_ = 0.005;
        }
        stan::test::expect_ad<true>(tols, f, phi_dbl[0], phi_dbl[1]);
      },
      theta_0);
}
