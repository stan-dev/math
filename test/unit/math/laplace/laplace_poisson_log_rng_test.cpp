#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

#include <gtest/gtest.h>

TEST_F(laplace_count_two_dim_diag_test, poisson_log_likelihood) {
  using stan::math::laplace_latent_poisson_log_rng;
  using stan::math::laplace_latent_tol_poisson_log_rng;
  using stan::math::multi_normal_rng;
  using stan::math::sqrt;
  using stan::math::square;

  // Compute exact mean and covariance.
  Eigen::VectorXd theta_root = stan::math::algebra_solver(
      stan::math::test::stationary_point(), theta_0, phi, d0, di0);
  Eigen::MatrixXd K_laplace
      = stan::math::test::laplace_covariance(theta_root, phi);

  boost::random::mt19937 rng;
  rng.seed(1954);
  Eigen::MatrixXd theta_pred = laplace_latent_poisson_log_rng(
      y, y_index, 0, stan::math::test::diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);

  // double tol = 1e-3;
  EXPECT_NEAR(theta_benchmark(0), theta_pred(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred(1), tol);

  // int n_sim = 5e5;
  Eigen::VectorXd theta_dim0(n_sim);
  Eigen::VectorXd theta_dim1(n_sim);
  for (int i = 0; i < n_sim; i++) {
    rng.seed(2025 + i);
    Eigen::MatrixXd theta_pred = laplace_latent_poisson_log_rng(
        y, y_index, 0, stan::math::test::diagonal_kernel_functor{},
        std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);

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
  using stan::math::laplace_latent_poisson_log_rng;
  using stan::math::log;
  using stan::math::multi_normal_rng;
  using stan::math::sqrt;
  using stan::math::square;

  rng.seed(1954);
  Eigen::MatrixXd theta_pred_exp = laplace_latent_poisson_log_rng(
      y, y_index, log(ye), stan::math::test::diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);

  EXPECT_NEAR(theta_benchmark(0), theta_pred_exp(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred_exp(1), tol);

  Eigen::VectorXd theta_dim0(n_sim);
  Eigen::VectorXd theta_dim1(n_sim);
  for (int i = 0; i < n_sim; i++) {
    rng.seed(2025 + i);
    Eigen::MatrixXd theta_pred = laplace_latent_poisson_log_rng(
        y, y_index, log(ye), stan::math::test::diagonal_kernel_functor{},
        std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);

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
