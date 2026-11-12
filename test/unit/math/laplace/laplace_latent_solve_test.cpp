#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/util.hpp>

#include <boost/random/mersenne_twister.hpp>

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
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
  auto [mean_est, chol_est]
      = laplace_latent_solve(poisson_log_likelihood{}, std::forward_as_tuple(y),
                             1, stan::math::test::diagonal_kernel_functor{},
                             std::forward_as_tuple(phi(0), phi(1)), nullptr);
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

TEST_F(laplace_count_two_dim_diag_test, latent_solve_tol_mean_and_cov) {
  using stan::math::laplace_latent_solve_tol;
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;
  auto [mean_est, chol_est] = laplace_latent_solve_tol(
      poisson_log_likelihood{}, std::forward_as_tuple(y), hessian_block_size,
      stan::math::test::diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                      max_steps_line_search, true),
      nullptr);
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
                     nullptr);
               }),
               std::domain_error);
}

namespace {
// y ~ normal(theta, s): the log likelihood is quadratic in theta, so the
// Laplace approximation is exact,
//   Sigma = (K^-1 + I / s^2)^-1,  mean = Sigma * y / s^2.
struct normal_likelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, const Eigen::VectorXd& y, double s,
                  std::ostream* pstream) const {
    return stan::math::normal_lpdf(y, theta, s);
  }
};

// y ~ cauchy(theta, s): not log-concave.  The negative Hessian of the log
// likelihood is diagonal with
//   W_j = (2 / s^2) (1 - r_j^2) / (1 + r_j^2)^2,  r_j = (y_j - theta_j) / s,
// which is negative for |r_j| > 1.
struct cauchy_likelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, const Eigen::VectorXd& y, double s,
                  std::ostream* pstream) const {
    return stan::math::cauchy_lpdf(y, theta, s);
  }
};

struct fixed_covariance {
  Eigen::MatrixXd operator()(const Eigen::MatrixXd& K,
                             std::ostream* pstream) const {
    return K;
  }
};

Eigen::MatrixXd prior_covariance() {
  Eigen::MatrixXd K(2, 2);
  K << 1.0, 0.5, 0.5, 2.0;
  return K;
}

// Laplace covariance (K^-1 + diag(W))^-1
Eigen::MatrixXd laplace_covariance(const Eigen::MatrixXd& K,
                                   const Eigen::VectorXd& W) {
  Eigen::MatrixXd precision = K.inverse();
  precision.diagonal() += W;
  return precision.inverse();
}

Eigen::VectorXd cauchy_neg_hessian(const Eigen::VectorXd& y,
                                   const Eigen::VectorXd& theta, double s) {
  Eigen::ArrayXd r2 = ((y - theta) / s).array().square();
  return ((2.0 / (s * s)) * (1.0 - r2) / (1.0 + r2).square()).matrix();
}

std::vector<Eigen::VectorXd> cauchy_data() {
  std::vector<Eigen::VectorXd> ys(2, Eigen::VectorXd(2));
  ys[0] << 2.5, -2.0;
  ys[1] << 4.0, -3.0;
  return ys;
}

template <typename LL>
auto latent_solve_tol(const LL& ll, const Eigen::VectorXd& y, double s,
                      int hessian_block_size, int solver,
                      bool allow_fallthrough) {
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(2);
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  constexpr int max_steps_line_search = 1000;
  return stan::math::laplace_latent_solve_tol(
      ll, std::forward_as_tuple(y, s), hessian_block_size, fixed_covariance{},
      std::make_tuple(prior_covariance()),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                      max_steps_line_search, allow_fallthrough),
      nullptr);
}
}  // namespace

// All solvers must return the same Laplace covariance.  Solver 2 used to
// return K - K W B^-1 W K instead of K_root B^-1 K_root^T.
TEST(laplace_latent_solve, all_solvers_match_closed_form_normal) {
  const Eigen::MatrixXd K = prior_covariance();
  Eigen::VectorXd y(2);
  y << 0.3, -0.2;
  const double s = 2.0;
  const Eigen::MatrixXd Sigma
      = laplace_covariance(K, Eigen::VectorXd::Constant(2, 1.0 / (s * s)));
  const Eigen::VectorXd mean = Sigma * y / (s * s);
  for (int hessian_block_size : {1, 2}) {
    for (int solver : {1, 2, 3}) {
      SCOPED_TRACE("hessian_block_size = " + std::to_string(hessian_block_size)
                   + ", solver = " + std::to_string(solver));
      auto [mean_est, chol_est] = latent_solve_tol(
          normal_likelihood{}, y, s, hessian_block_size, solver, false);
      EXPECT_MATRIX_NEAR(mean, mean_est, 1e-8);
      EXPECT_MATRIX_NEAR(Sigma, (chol_est * chol_est.transpose()).eval(), 1e-8);
    }
  }
}

// With allow_fallthrough, solver 1 fails on a likelihood that is not
// log-concave, and the solve continues with solver 2 or 3.  The covariance
// must come from the solver that converged.  The solver-1 formula used to be
// applied: a wrong covariance after a fallthrough to solver 2, and a failed
// Eigen assertion (L is 0 x 0) after a fallthrough to solver 3.
TEST(laplace_latent_solve, fallthrough_matches_closed_form_cauchy) {
  const Eigen::MatrixXd K = prior_covariance();
  const double s = 0.5;
  for (const auto& y : cauchy_data()) {
    // solver 3 handles an indefinite W without a fallthrough
    auto [mean_ref, chol_ref]
        = latent_solve_tol(cauchy_likelihood{}, y, s, 2, 3, false);
    const Eigen::VectorXd W = cauchy_neg_hessian(y, mean_ref, s);
    ASSERT_LT(W.minCoeff(), 0.0);  // not log-concave at the mode
    const Eigen::MatrixXd Sigma = laplace_covariance(K, W);
    EXPECT_MATRIX_NEAR(Sigma, (chol_ref * chol_ref.transpose()).eval(), 1e-6);
    for (int solver : {1, 2}) {
      SCOPED_TRACE("y = (" + std::to_string(y(0)) + ", " + std::to_string(y(1))
                   + "), solver = " + std::to_string(solver));
      auto [mean_est, chol_est]
          = latent_solve_tol(cauchy_likelihood{}, y, s, 2, solver, true);
      EXPECT_MATRIX_NEAR(mean_ref, mean_est, 1e-6);
      EXPECT_MATRIX_NEAR(Sigma, (chol_est * chol_est.transpose()).eval(), 1e-6);
    }
  }
}

// The functions without control arguments use solver 1 with
// allow_fallthrough, so a likelihood that is not log-concave takes the same
// path.  The default tolerance (1.49e-8 on the change of the objective)
// gives a mode error of order sqrt(1.49e-8), hence the test tolerance 1e-3.
TEST(laplace_latent_solve, default_options_cauchy) {
  const Eigen::MatrixXd K = prior_covariance();
  const double s = 0.5;
  for (const auto& y : cauchy_data()) {
    auto [mean_ref, chol_ref]
        = latent_solve_tol(cauchy_likelihood{}, y, s, 2, 3, false);
    const Eigen::MatrixXd Sigma
        = laplace_covariance(K, cauchy_neg_hessian(y, mean_ref, s));
    auto [mean_est, chol_est] = stan::math::laplace_latent_solve(
        cauchy_likelihood{}, std::forward_as_tuple(y, s), 2, fixed_covariance{},
        std::make_tuple(K), nullptr);
    EXPECT_MATRIX_NEAR(mean_ref, mean_est, 1e-3);
    EXPECT_MATRIX_NEAR(Sigma, (chol_est * chol_est.transpose()).eval(), 1e-3);
  }
}

// laplace_latent_tol_rng() shares the covariance computation: with the same
// seed, its draw equals multi_normal_rng() from the exact mean and
// covariance.
TEST(laplace_latent_tol_rng, solver_2_draw_uses_laplace_covariance) {
  const Eigen::MatrixXd K = prior_covariance();
  Eigen::VectorXd y(2);
  y << 0.3, -0.2;
  const double s = 2.0;
  const Eigen::MatrixXd Sigma
      = laplace_covariance(K, Eigen::VectorXd::Constant(2, 1.0 / (s * s)));
  const Eigen::VectorXd mean = Sigma * y / (s * s);
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(2);
  boost::random::mt19937 rng_laplace(1234);
  boost::random::mt19937 rng_ref(1234);
  Eigen::VectorXd draw = stan::math::laplace_latent_tol_rng(
      normal_likelihood{}, std::forward_as_tuple(y, s), 2, fixed_covariance{},
      std::make_tuple(K), std::make_tuple(theta_0, 1e-12, 1000, 2, 1000, false),
      rng_laplace, nullptr);
  Eigen::VectorXd draw_ref = stan::math::multi_normal_rng(mean, Sigma, rng_ref);
  EXPECT_MATRIX_NEAR(draw_ref, draw, 1e-8);
}
