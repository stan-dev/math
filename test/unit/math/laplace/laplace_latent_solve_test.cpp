#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>

#include <boost/random/mersenne_twister.hpp>

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

namespace {
struct poisson_log_likelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, const std::vector<int>& y,
                  std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(y, theta);
  }
};
}  // namespace

TEST_F(laplace_count_two_dim_diag_test, latent_solve_mean_and_cov) {
  using stan::math::laplace_latent_solve;
  auto [mean_est, chol_est] = laplace_latent_solve(
      poisson_log_likelihood{}, std::forward_as_tuple(y), 1,
      stan::math::test::diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);
  constexpr double tol = 1e-6;
  EXPECT_EQ(2, mean_est.size());
  EXPECT_NEAR(theta_root(0), mean_est(0), tol);
  EXPECT_NEAR(theta_root(1), mean_est(1), tol);
  EXPECT_NEAR(0.0, chol_est(0, 1), 1e-12);  // check lower triangular matrix
  Eigen::MatrixXd Sigma_est = chol_est * chol_est.transpose();
  EXPECT_NEAR(K_laplace(0, 0), Sigma_est(0, 0), tol);
  EXPECT_NEAR(K_laplace(1, 1), Sigma_est(1, 1), tol);
  EXPECT_NEAR(K_laplace(0, 1), Sigma_est(0, 1), tol);
  EXPECT_NEAR(K_laplace(1, 0), Sigma_est(1, 0), tol);
}

TEST_F(laplace_count_two_dim_diag_test, latent_tol_solve_mean_and_cov) {
  using stan::math::laplace_latent_tol_solve;
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  auto [mean_est, chol_est] = laplace_latent_tol_solve(
      poisson_log_likelihood{}, std::forward_as_tuple(y), hessian_block_size,
      stan::math::test::diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                      max_steps_line_search, true),
      rng, nullptr);
  constexpr double tol = 1e-6;
  EXPECT_EQ(2, mean_est.size());
  EXPECT_NEAR(theta_root(0), mean_est(0), tol);
  EXPECT_NEAR(theta_root(1), mean_est(1), tol);
  EXPECT_NEAR(0.0, chol_est(0, 1), 1e-12);  // check lower triangular matrix
  Eigen::MatrixXd Sigma_est = chol_est * chol_est.transpose();
  EXPECT_NEAR(K_laplace(0, 0), Sigma_est(0, 0), tol);
  EXPECT_NEAR(K_laplace(1, 1), Sigma_est(1, 1), tol);
  EXPECT_NEAR(K_laplace(0, 1), Sigma_est(0, 1), tol);
  EXPECT_NEAR(K_laplace(1, 0), Sigma_est(1, 0), tol);
}

TEST_F(laplace_count_two_dim_diag_test,
       latent_solve_singular_covariance_throws) {
  using stan::math::laplace_latent_solve;
  EXPECT_THROW(({
                 laplace_latent_solve(
                     poisson_log_likelihood{}, std::forward_as_tuple(y), 1,
                     stan::math::test::diagonal_kernel_functor{},
                     std::forward_as_tuple(0.0, phi(1)),  // singular covariance
                     rng, nullptr);
               }),
               std::domain_error);
}
