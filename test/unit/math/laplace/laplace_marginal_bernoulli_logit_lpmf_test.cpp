#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
// #include <stan/math/laplace/laplace_likelihood_bernoulli_logit.hpp>
// #include <stan/math/laplace/laplace_marginal_bernoulli_logit_lpmf.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>

#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

// TEST(laplace, likelihood_differentiation) {
//   using stan::math::diff_bernoulli_logit;
//   using stan::math::var;
//
//   double test_tolerance = 2e-4;
//
//   Eigen::VectorXd theta(2);
//   theta << -2.45809, -3.6127;
//   Eigen::VectorXd y(2), n_samples(2);
//   y << 1, 0;
//   n_samples << 1, 1;
//   Eigen::Matrix<var, Eigen::Dynamic, 1> theta_v = theta;
//   Eigen::VectorXd eta_dummy;
//
//   diff_bernoulli_logit diff_functor(n_samples, y);
//   double log_density = diff_functor.log_likelihood(theta, eta_dummy);
//   Eigen::VectorXd gradient;
//   Eigen::SparseMatrix<double> hessian;
//   diff_functor.diff(theta, eta_dummy, gradient, hessian);
//   Eigen::VectorXd third_tensor = diff_functor.third_diff(theta, eta_dummy);
//
//   EXPECT_NEAR(-2.566843, log_density, test_tolerance);
//
//   // finite diff calculations for first-order derivatives
//   double diff = 1e-12;
//   Eigen::VectorXd theta_1u = theta;
//   Eigen::VectorXd theta_1l = theta;
//   Eigen::VectorXd theta_2u = theta;
//   Eigen::VectorXd theta_2l = theta;
//   theta_1u(0) = theta(0) + diff;
//   theta_1l(0) = theta(0) - diff;
//   theta_2u(1) = theta(1) + diff;
//   theta_2l(1) = theta(1) - diff;
//   double diff_1 = (diff_functor.log_likelihood(theta_1u, eta_dummy)
//                    - diff_functor.log_likelihood(theta_1l, eta_dummy))
//                   / (2 * diff);
//   double diff_2 = (diff_functor.log_likelihood(theta_2u, eta_dummy)
//                    - diff_functor.log_likelihood(theta_2l, eta_dummy))
//                   / (2 * diff);
//
//   EXPECT_NEAR(diff_1, gradient(0), test_tolerance);
//   EXPECT_NEAR(diff_2, gradient(1), test_tolerance);
//
//   // finite diff calculation for second-order derivatives
//   Eigen::VectorXd gradient_1u, gradient_1l, gradient_2u, gradient_2l;
//   Eigen::SparseMatrix<double> hessian_1u, hessian_1l, hessian_2u, hessian_2l;
//   diff_functor.diff(theta_1u, eta_dummy, gradient_1u, hessian_1u);
//   diff_functor.diff(theta_1l, eta_dummy, gradient_1l, hessian_1l);
//   diff_functor.diff(theta_2u, eta_dummy, gradient_2u, hessian_2u);
//   diff_functor.diff(theta_2l, eta_dummy, gradient_2l, hessian_2l);
//
//   double diff_grad_1 = (gradient_1u(0) - gradient_1l(0)) / (2 * diff);
//   double diff_grad_2 = (gradient_2u(1) - gradient_2l(1)) / (2 * diff);
//
//   EXPECT_NEAR(diff_grad_1, hessian.coeff(0, 0), test_tolerance);
//   EXPECT_NEAR(diff_grad_2, hessian.coeff(1, 1), test_tolerance);
//
//   // finite diff calculation for third-order derivatives
//   double diff_hess_1
//       = (hessian_1u.coeff(0, 0) - hessian_1l.coeff(0, 0)) / (2 * diff);
//   double diff_hess_2
//       = (hessian_2u.coeff(1, 1) - hessian_2l.coeff(1, 1)) / (2 * diff);
//
//   EXPECT_NEAR(diff_hess_1, third_tensor(0), test_tolerance);
//   EXPECT_NEAR(diff_hess_2, third_tensor(1), test_tolerance);
// }

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

  Eigen::VectorXd theta_laplace, a, l_grad;
  Eigen::SparseMatrix<double> W_root;
  Eigen::MatrixXd L, covariance;
  std::vector<double> delta;
  std::vector<int> delta_int;
  int dim_phi = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
  phi << 1.6, 1;

  stan::math::test::sqr_exp_kernel_functor K;
  var target = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, K, phi, x,
                                                     delta, delta_int, theta_0);

  double tol = 8e-5;
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, value_of(target), tol);

  // Test with optional arguments.
  double tolerance = 1e-6;
  int max_num_steps = 100;
  target = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, K, phi, x,
                                                 delta, delta_int, theta_0,
                                                 0, tolerance, max_num_steps);
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

  double target_1u = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, K,
                                          phi_1u, x, delta, delta_int, theta_0),
         target_1l = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, K,
                                          phi_1l, x, delta, delta_int, theta_0),
         target_2u = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, K,
                                          phi_2u, x, delta, delta_int, theta_0),
         target_2l = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, K,
                                          phi_2l, x, delta, delta_int, theta_0);

  std::vector<double> g_finite(dim_phi);
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);

  tol = 1e-3;
  EXPECT_NEAR(g_finite[0], g[0], tol);
  EXPECT_NEAR(g_finite[1], g[1], tol);

}
