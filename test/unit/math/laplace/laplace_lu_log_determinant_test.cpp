#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/util.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>

namespace {
struct normal_likelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, const Eigen::VectorXd& y, double s,
                  std::ostream* pstream) const {
    return stan::math::normal_lpdf(y, theta, s);
  }
};

// theta(0) has two Cauchy observations y, which make log p(theta | y)
// bimodal in theta(0); theta(1) has one normal observation z.
struct cauchy_normal_likelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, const Eigen::VectorXd& y, double s,
                  double z, double sz, std::ostream* pstream) const {
    return stan::math::cauchy_lpdf(y, theta(0), s)
           + stan::math::normal_lpdf(z, theta(1), sz);
  }
};

struct fixed_covariance {
  Eigen::MatrixXd operator()(const Eigen::MatrixXd& K,
                             std::ostream* pstream) const {
    return K;
  }
};
}  // namespace

// Solver 3 factorizes B = I + K W with partial pivoting.  Here the rows are
// swapped and U has a negative diagonal entry although det(B) = 226 > 0.  The
// log determinant used to be sum(log(diag(U))), which is NaN.  For a normal
// likelihood the Laplace approximation is exact, so all solvers must return
// log N(y | 0, K + s^2 I).
TEST(laplace_marginal_tol, all_solvers_exact_with_negative_lu_pivot) {
  Eigen::MatrixXd K(2, 2);
  K << 1.0, 1.5, 1.5, 4.0;
  Eigen::VectorXd y(2);
  y << 0.3, -0.2;
  const double s = std::sqrt(0.1);
  const Eigen::MatrixXd B
      = Eigen::MatrixXd::Identity(2, 2) + K / (s * s);  // W = I / s^2
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(B);
  ASSERT_LT(lu.matrixLU().diagonal().minCoeff(), 0.0);
  ASSERT_GT(B.determinant(), 0.0);
  const double exact = stan::math::multi_normal_lpdf(
      y, Eigen::VectorXd::Zero(2),
      (K + s * s * Eigen::MatrixXd::Identity(2, 2)).eval());
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(2);
  for (int solver : {1, 2, 3}) {
    SCOPED_TRACE("solver = " + std::to_string(solver));
    const double lml = stan::math::laplace_marginal_tol<false>(
        normal_likelihood{}, std::forward_as_tuple(y, s), 2, fixed_covariance{},
        std::make_tuple(K),
        std::make_tuple(theta_0, 1e-12, 1000, solver, 1000, false), nullptr);
    EXPECT_NEAR(exact, lml, 1e-8);
  }
}

// Started near the saddle point (0, 0) of log p(theta | y), the Newton
// iteration of solver 3 converges to it.  There Sigma^-1 + W is not positive
// definite and the Laplace approximation is not defined.  This must be
// reported as such: laplace_marginal_tol() used to return NaN or a finite
// value, and laplace_latent_solve_tol() threw from cholesky_decompose().
// Started away from the saddle point, solver 3 finds a maximum.
TEST(laplace_marginal_tol, solver_3_saddle_point_throws) {
  Eigen::MatrixXd K(2, 2);
  K << 4.0, 1.0, 1.0, 1.0;
  Eigen::VectorXd y(2);
  y << -3.0, 3.0;
  const double s = 0.3;
  const double z = 0.0;
  const double sz = 0.1;
  const std::string reason = "det(I + Sigma * W) is not positive";
  Eigen::VectorXd theta_near(2);
  theta_near << 0.01, 3.0;
  auto ops_near = std::make_tuple(theta_near, 1e-10, 100, 3, 100, false);
  try {
    stan::math::laplace_marginal_tol<false>(
        cauchy_normal_likelihood{}, std::forward_as_tuple(y, s, z, sz), 2,
        fixed_covariance{}, std::make_tuple(K), ops_near, nullptr);
    FAIL() << "laplace_marginal_tol: expected std::domain_error";
  } catch (const std::domain_error& e) {
    EXPECT_NE(std::string(e.what()).find(reason), std::string::npos)
        << e.what();
  }
  try {
    stan::math::laplace_latent_solve_tol(
        cauchy_normal_likelihood{}, std::forward_as_tuple(y, s, z, sz), 2,
        fixed_covariance{}, std::make_tuple(K), ops_near, nullptr);
    FAIL() << "laplace_latent_solve_tol: expected std::domain_error";
  } catch (const std::domain_error& e) {
    EXPECT_NE(std::string(e.what()).find(reason), std::string::npos)
        << e.what();
  }
  Eigen::VectorXd theta_away(2);
  theta_away << 0.3, 3.0;
  auto [mean, chol] = stan::math::laplace_latent_solve_tol(
      cauchy_normal_likelihood{}, std::forward_as_tuple(y, s, z, sz), 2,
      fixed_covariance{}, std::make_tuple(K),
      std::make_tuple(theta_away, 1e-10, 100, 3, 100, false), nullptr);
  EXPECT_NEAR(2.9385, mean(0), 1e-4);
  EXPECT_NEAR(0.0097, mean(1), 1e-4);
}
