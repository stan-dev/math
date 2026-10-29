#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <utility>
#include <vector>

TEST_F(AgradRev, mathMixScalFun_neg_binomial_2_log_glm_lpmf) {
  auto f = [](const auto y, const auto& x) {
    return [=](const auto& alpha, const auto& beta, const auto& phi) {
      return stan::math::neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi);
    };
  };
  auto f2 = [](const auto y, const auto& phi) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi);
    };
  };

  std::vector<int> y{0, 1};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);
  double phi = 1.5;

  stan::test::expect_ad(f(y[0], x), alpha, beta, phi);
  stan::test::expect_ad(f(y[0], x), alpha[0], beta, phi);
  stan::test::expect_ad(f(y[0], x_rowvec), alpha, beta, phi);
  stan::test::expect_ad(f(y[0], x_rowvec), alpha[0], beta, phi);

  stan::test::expect_ad(f2(y, phi), x, alpha, beta);
  stan::test::expect_ad(f2(y, phi), x, alpha[0], beta);
  stan::test::expect_ad(f2(y, phi), x_rowvec, alpha, beta);
  stan::test::expect_ad(f2(y, phi), x_rowvec, alpha[0], beta);
}

// x is a column vector: a design matrix with a single attribute
TEST_F(AgradRev, mathMixScalFun_neg_binomial_2_log_glm_lpmf_x_column_vector) {
  std::vector<int> y{0, 1};
  auto f = [y](const auto& x, const auto& alpha, const auto& beta) {
    return stan::math::neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, 1.5);
  };

  Eigen::VectorXd x = Eigen::VectorXd::Random(2);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(1);

  stan::test::expect_ad(f, x, alpha, beta);
}

// The second derivative in the linear predictor theta at large theta,
// -(y + phi) phi mu / (mu + phi)^2 with mu = exp(theta). A quotient-rule form
// of the theta partial cancels there. References from mpmath at 100 digits,
// the closed form checked against mp.diff.
TEST_F(AgradRev,
       mathMixScalFun_neg_binomial_2_log_glm_lpmf_large_theta_hessian) {
  using stan::math::fvar;
  using stan::math::var;
  const std::vector<int> y{1};
  Eigen::MatrixXd x(1, 1);
  x << 1.0;
  Eigen::VectorXd beta(1);
  beta << 0.0;
  const double phi = 0.07;
  const std::vector<std::pair<double, double>> cases
      = {{35.0, -4.7225324533500949e-17}, {27.0, -1.4077670835874031e-13}};
  for (const auto& c : cases) {
    fvar<var> alpha(var(c.first), 1.0);
    fvar<var> lp
        = stan::math::neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi);
    lp.d_.grad();
    EXPECT_NEAR(alpha.val_.adj(), c.second, 1e-10 * std::fabs(c.second))
        << "theta = " << c.first;
    stan::math::recover_memory();
  }
}
