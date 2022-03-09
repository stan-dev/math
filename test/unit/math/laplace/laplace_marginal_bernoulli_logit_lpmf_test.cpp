#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>

#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace_marginal_bernoulli_logit_lpmf, phi_dim500) {
  using stan::math::diff_bernoulli_logit;
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;
  using stan::math::to_vector;
  using stan::math::var;

  int dim_theta = 500;
  int n_observations = 500;
  std::string data_directory = "test/unit/math/laplace/aki_synth_data/";
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
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
  phi << 1.6, 1;

  stan::math::test::sqr_exp_kernel_functor kernel_functor;
  var target = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                                   nullptr, x, phi(0), phi(1));

  double tol = 8e-5;
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, value_of(target), tol);

  // Test with optional arguments.
  double tolerance = 1e-6;
  int max_num_steps = 100;
  target = laplace_marginal_tol_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                                 0, tolerance, max_num_steps, 0, 1, 0, x, phi(0), phi(1));
  EXPECT_NEAR(-195.368, value_of(target), tol);

  std::vector<double> g;
  std::vector<stan::math::var> parm_vec{phi(0), phi(1)};
  target.grad(parm_vec, g);

  // finite diff benchmark
  double diff = 1e-7;
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_1l = phi_dbl, phi_1u = phi_dbl, phi_2l = phi_dbl,
                  phi_2u = phi_dbl;
  phi_1l(0) -= diff;
  phi_1u(0) += diff;
  phi_2l(1) -= diff;
  phi_2u(1) += diff;

  double target_1u = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                                   nullptr, x, phi_1u(0), phi_1u(1));
  double target_1l = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                                   nullptr, x, phi_1l(0), phi_1l(1));
  double target_2u = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                                   nullptr, x, phi_2u(0), phi_2u(1));
  double target_2l = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                                   nullptr, x, phi_2l(0), phi_2l(1));

  std::vector<double> g_finite(dim_phi);
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);

  tol = 1e-3;
  EXPECT_NEAR(g_finite[0], g[0], tol);
  EXPECT_NEAR(g_finite[1], g[1], tol);
}
