#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/mix/laplace/laplace_utility.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST_F(laplace_count_two_dim_diag_test, poisson_log_likelihood) {
  using stan::math::laplace_marginal_poisson_log_rng;
  using stan::math::laplace_marginal_tol_poisson_log_rng;
  using stan::math::multi_normal_rng;
  using stan::math::sqrt;
  using stan::math::square;

  // Compute exact mean and covariance.
  Eigen::VectorXd theta_root
      = stan::math::algebra_solver(stationary_point(), theta_0, phi, d0, di0);
  Eigen::MatrixXd K_laplace = laplace_covariance(theta_root, phi);

  boost::random::mt19937 rng;
  rng.seed(1954);
  Eigen::MatrixXd theta_pred = laplace_marginal_poisson_log_rng(
      y, y_index, theta_0, diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)), std::make_tuple(),
      std::make_tuple(), rng, nullptr);

  // double tol = 1e-3;
  EXPECT_NEAR(theta_benchmark(0), theta_pred(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred(1), tol);

  // int n_sim = 5e5;
  Eigen::VectorXd theta_dim0(n_sim);
  Eigen::VectorXd theta_dim1(n_sim);
  for (int i = 0; i < n_sim; i++) {
    rng.seed(2025 + i);
    Eigen::MatrixXd theta_pred = laplace_marginal_poisson_log_rng(
        y, y_index, theta_0, diagonal_kernel_functor{},
        std::forward_as_tuple(phi(0), phi(1)), std::make_tuple(),
        std::make_tuple(), rng, nullptr);

    theta_dim0(i) = theta_pred(0);
    theta_dim1(i) = theta_pred(1);
  }

  Eigen::MatrixXd K_sample(2, 2);
  K_sample(0, 0)
      = theta_dim0.array().square().mean() - square(theta_dim0.mean());
  K_sample(1, 1)
      = theta_dim1.array().square().mean() - square(theta_dim1.mean());
  K_sample(0, 1) = theta_dim0.cwiseProduct(theta_dim1).mean()
                   - theta_dim0.mean() * theta_dim1.mean();
  K_sample(1, 0) = K_sample(0, 1);

  // Check answers are within three std of the true answer.
  EXPECT_NEAR(theta_root(0), theta_dim0.mean(),
              3 * sqrt(K_laplace(0, 0) / n_sim));
  EXPECT_NEAR(theta_root(1), theta_dim1.mean(),
              3 * sqrt(K_laplace(1, 1) / n_sim));

  // Check sample covariance
  EXPECT_NEAR(K_laplace(0, 0), K_sample(0, 0), 5e-3);
  EXPECT_NEAR(K_laplace(1, 1), K_sample(1, 1), 6e-3);
  EXPECT_NEAR(K_laplace(0, 1), K_sample(0, 1), 1e-3);
}

TEST_F(laplace_count_two_dim_diag_test, poisson_log_exp_likelihood) {
  using stan::math::laplace_marginal_poisson_2_log_rng;
  using stan::math::multi_normal_rng;
  using stan::math::sqrt;
  using stan::math::square;

  rng.seed(1954);
  Eigen::MatrixXd theta_pred_exp = laplace_marginal_poisson_2_log_rng(
      y, y_index, ye, theta_0, diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)), std::make_tuple(),
      std::make_tuple(), rng, nullptr);

  EXPECT_NEAR(theta_benchmark(0), theta_pred_exp(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred_exp(1), tol);

  Eigen::VectorXd theta_dim0(n_sim);
  Eigen::VectorXd theta_dim1(n_sim);
  for (int i = 0; i < n_sim; i++) {
    rng.seed(2025 + i);
    Eigen::MatrixXd theta_pred = laplace_marginal_poisson_2_log_rng(
        y, y_index, ye, theta_0, diagonal_kernel_functor{},
        std::forward_as_tuple(phi(0), phi(1)), std::make_tuple(),
        std::make_tuple(), rng, nullptr);

    theta_dim0(i) = theta_pred(0);
    theta_dim1(i) = theta_pred(1);
  }

  Eigen::MatrixXd K_sample(2, 2);
  K_sample(0, 0)
      = theta_dim0.array().square().mean() - square(theta_dim0.mean());
  K_sample(1, 1)
      = theta_dim1.array().square().mean() - square(theta_dim1.mean());
  K_sample(0, 1) = theta_dim0.cwiseProduct(theta_dim1).mean()
                   - theta_dim0.mean() * theta_dim1.mean();
  K_sample(1, 0) = K_sample(0, 1);

  // Check answers are within three std of the true answer.
  EXPECT_NEAR(theta_root(0), theta_dim0.mean(),
              3 * sqrt(K_laplace(0, 0) / n_sim));
  EXPECT_NEAR(theta_root(1), theta_dim1.mean(),
              3 * sqrt(K_laplace(1, 1) / n_sim));

  // Check sample covariance
  EXPECT_NEAR(K_laplace(0, 0), K_sample(0, 0), 5e-3);
  EXPECT_NEAR(K_laplace(1, 1), K_sample(1, 1), 6e-3);
  EXPECT_NEAR(K_laplace(0, 1), K_sample(0, 1), 1e-3);
}

// Keep this or delete this?
// @charlesm93: no need to include this for now but there are some
// interesting ideas here.
/*
TEST(laplace, poisson_basic_rng) {
  using stan::math::algebra_solver;
  using stan::math::diag_matrix;
  using stan::math::diag_pre_multiply;
  using stan::math::diff_poisson_log;
  using stan::math::inv;
  using stan::math::laplace_poisson_log_rng;
  using stan::math::laplace_rng;
  using stan::math::mdivide_left_tri;
  using stan::math::square;
  using stan::math::to_vector;
  using stan::math::value_of;

  Eigen::VectorXd theta_0(2);
  theta_0 << 1, 1;
  Eigen::VectorXd sigma(2);
  sigma << 3, 2;
  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};
  stan::math::poisson_log_likelihood laplace_likelihood{};
//  diff_poisson_log laplace_likelihood(to_vector(n_samples), to_vector(sums));
  std::vector<double> d0;
  std::vector<int> di0;

  // Method 1: brute force and straightforward
  Eigen::VectorXd theta_root
      = algebra_solver(stationary_point(), theta_0, sigma, d0, di0);

  Eigen::VectorXd gradient;
  Eigen::SparseMatrix<double> W_sparse;
  Eigen::VectorXd eta_dummy;
  std::tie(gradient, eta_dummy, W_sparse) =
stan::math::laplace_likelihood::internal::diff(laplace_likelihood, theta_root);
  Eigen::MatrixXd W = -W_sparse;
  diagonal_kernel_functor covariance_function;
  std::vector<Eigen::VectorXd> x_dummy;
  Eigen::MatrixXd K = covariance_function(sigma, x_dummy, d0, di0, 0);

  std::cout << "K (brute force): " << std::endl
            << (K.inverse() + W).inverse() << std::endl
            << std::endl;

  // Method 2: Vectorized R&W method
  constexpr double tolerance = 1e-8;
  constexpr int max_num_steps = 1000;
  Eigen::MatrixXd K_root;
  // First find the mode using the custom Newton step
  Eigen::MatrixXd covariance;
  Eigen::VectorXd theta;
  Eigen::SparseMatrix<double> W_r;
  Eigen::MatrixXd L;
  Eigen::VectorXd theta0_val = value_of(theta_0);
  {
    Eigen::VectorXd a;
    Eigen::VectorXd l_grad;
    Eigen::PartialPivLU<Eigen::MatrixXd> LU_dummy;
    double marginal_density = laplace_marginal_density(
        laplace_likelihood, std::forward_as_tuple(sums, n_samples),
covariance_function, std::forward_as_tuple(sigma, d0, di0), covariance, theta,
W_r, L, a, l_grad, LU_dummy, K_root, theta0_val, 0, tolerance, max_num_steps);

    std::cout << "theta (mode) = " << theta.transpose() << std::endl;
  }

  Eigen::MatrixXd V;
  Eigen::VectorXd W_root(theta.size());
  for (int i = 0; i < theta.size(); i++)
    W_root(i) = W_r.coeff(i, i);
  V = mdivide_left_tri<Eigen::Lower>(L, diag_pre_multiply(W_root, covariance));
  std::cout << "K (method 1): " << std::endl
            << covariance - V.transpose() * V << std::endl
            << std::endl;

  // Method 3: Modified R&W method
  Eigen::VectorXd W_root_inv = inv(W_root);
  Eigen::MatrixXd V_dec
      = mdivide_left_tri<Eigen::Lower>(L, diag_matrix(W_root_inv));
  std::cout << "K (method 2): " << std::endl
            << -V_dec.transpose() * V_dec + diag_matrix(square(W_root_inv))
            << std::endl
            << std::endl;

  // Call to rng function
  boost::random::mt19937 rng;
  Eigen::MatrixXd theta_pred
      = laplace_poisson_log_rng(sums, n_samples, covariance_function, sigma,
                                x_dummy, d0, di0, theta_0, rng);
}
*/
