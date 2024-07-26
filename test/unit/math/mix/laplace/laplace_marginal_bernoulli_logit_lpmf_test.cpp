#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/mix/laplace/laplace_utility.hpp>
#include <test/unit/math/mix/laplace/aki_synth_data/x1.hpp>

#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace_marginal_bernoulli_logit_lpmf, phi_dim500) {
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;
  using stan::math::laplace_marginal_tol_bernoulli_logit_lpmf;
  using stan::math::to_vector;
  using stan::math::var;

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
