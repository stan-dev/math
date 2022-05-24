#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>

#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>


TEST_F(laplace_bernoulli_dim500_test, dense_covariance) {
  using stan::math::diff_bernoulli_logit;
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;
  using stan::math::to_vector;
  using stan::math::var;

  stan::math::test::sqr_exp_kernel_functor kernel_functor;
  var target = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                                   nullptr, x, phi(0), phi(1));
  double tol = 8e-5;
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, value_of(target), tol);

  // Test with optional arguments.
  double tolerance = 1e-6;
  int max_num_steps = 100;
  target = laplace_marginal_tol_bernoulli_logit_lpmf(y, n_samples,
                    tolerance, max_num_steps, 1, 0, 1, 0,
                    theta_0, kernel_functor, 0, x, phi(0), phi(1));
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

TEST_F(laplace_bernoulli_dim500_test, diagonal_covariance) {
    using stan::math::laplace_marginal_bernoulli_logit_lpmf;
    using stan::math::laplace_marginal_tol_bernoulli_logit_lpmf;
    using stan::math::var;

    stan::math::test::sqr_exp_kernel_diag_functor kernel_functor;

    var target_dense =
      laplace_marginal_bernoulli_logit_lpmf(y, n_samples, theta_0, kernel_functor,
                                            nullptr, x, phi(0), phi(1),
                                            n_observations);

    double tolerance = 1e-6;
    int max_num_steps = 100;
    int hessian_block_size = 1;
    int covariance_block_size = 1;
    int solver = 1;
    int max_steps_line_search = 0;

    var target =
      laplace_marginal_tol_bernoulli_logit_lpmf(y, n_samples, tolerance,
                                            max_num_steps,
                                            hessian_block_size,
                                            covariance_block_size,
                                            solver,
                                            max_steps_line_search,
                                            theta_0,
                                            kernel_functor, nullptr,
                                            x, phi(0), phi(1),
                                            n_observations);

   double tol = 1e-6;
   EXPECT_NEAR(value_of(target_dense), value_of(target), tol);

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

   double target_1u = laplace_marginal_tol_bernoulli_logit_lpmf(y, n_samples,
    tolerance,
    max_num_steps,
    hessian_block_size,
    covariance_block_size,
    solver,
    max_steps_line_search,
    theta_0, kernel_functor,
    nullptr, x, phi_1u(0), phi_1u(1), n_observations);
   double target_1l = laplace_marginal_tol_bernoulli_logit_lpmf(y, n_samples,
     tolerance,
     max_num_steps,
     hessian_block_size,
     covariance_block_size,
     solver,
     max_steps_line_search,
     theta_0, kernel_functor,
     nullptr, x, phi_1l(0), phi_1l(1), n_observations);
   double target_2u = laplace_marginal_tol_bernoulli_logit_lpmf(y, n_samples,
     tolerance,
     max_num_steps,
     hessian_block_size,
     covariance_block_size,
     solver,
     max_steps_line_search,
     theta_0, kernel_functor,
     nullptr, x, phi_2u(0), phi_2u(1), n_observations);
   double target_2l = laplace_marginal_tol_bernoulli_logit_lpmf(y, n_samples,
     tolerance,
     max_num_steps,
     hessian_block_size,
     covariance_block_size,
     solver,
     max_steps_line_search,
     theta_0, kernel_functor,
     nullptr, x, phi_2l(0), phi_2l(1), n_observations);

   std::vector<double> g_finite(dim_phi);
   g_finite[0] = (target_1u - target_1l) / (2 * diff);
   g_finite[1] = (target_2u - target_2l) / (2 * diff);

   tol = 5e-3;
   EXPECT_NEAR(g_finite[0], g[0], tol);
   EXPECT_NEAR(g_finite[1], g[1], tol);

}
