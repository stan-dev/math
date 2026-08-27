#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

namespace generalized_normal_lpdf_test {

// expect_ad takes at most three autodiff arguments, so each factory holds
// one of (y, mu, alpha, beta) primitive and rotates the other three.
auto f_y_mu_alpha(double beta) {
  return [beta](const auto& y, const auto& mu, const auto& alpha) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, beta);
  };
}

auto f_y_mu_beta(double alpha) {
  return [alpha](const auto& y, const auto& mu, const auto& beta) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, beta);
  };
}

auto f_y_alpha_beta(double mu) {
  return [mu](const auto& y, const auto& alpha, const auto& beta) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, beta);
  };
}

auto f_mu_alpha_beta(double y) {
  return [y](const auto& mu, const auto& alpha, const auto& beta) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, beta);
  };
}

// Every argument is under autodiff in at least one of the four rotations.
void expect_all_rotations(double y, double mu, double alpha, double beta) {
  stan::test::expect_ad(f_y_mu_alpha(beta), y, mu, alpha);
  stan::test::expect_ad(f_y_mu_beta(alpha), y, mu, beta);
  stan::test::expect_ad(f_y_alpha_beta(mu), y, alpha, beta);
  stan::test::expect_ad(f_mu_alpha_beta(y), mu, alpha, beta);
}

}  // namespace generalized_normal_lpdf_test

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_interior) {
  using generalized_normal_lpdf_test::expect_all_rotations;

  expect_all_rotations(1.3, 0.4, 1.7, 2.0);
  expect_all_rotations(-2.5, 0.8, 0.6, 3.0);
  expect_all_rotations(0.0, -1.1, 2.3, 1.4);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_beta_boundaries) {
  using generalized_normal_lpdf_test::expect_all_rotations;

  // beta == 1 is the double exponential, beta == 2 the normal.
  expect_all_rotations(1.7, 0.3, 1.2, 1.0);
  expect_all_rotations(1.7, 0.3, 1.2, 2.0);

  // beta < 1 gives a cusp at y == mu; beta large approaches the uniform.
  expect_all_rotations(1.7, 0.3, 1.2, 0.5);
  expect_all_rotations(0.4, 0.3, 1.2, 8.0);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_alpha_boundaries) {
  using generalized_normal_lpdf_test::expect_all_rotations;

  expect_all_rotations(0.9, 0.5, 0.05, 2.0);
  expect_all_rotations(0.9, 0.5, 40.0, 2.0);
  expect_all_rotations(0.9, 0.5, 40.0, 0.7);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_tails) {
  using generalized_normal_lpdf_test::expect_all_rotations;

  // Large (|y - mu| / alpha) ^ beta, where pow() dominates the value.
  expect_all_rotations(12.0, 0.0, 1.0, 2.0);
  expect_all_rotations(-9.5, 1.5, 0.8, 1.0);
  expect_all_rotations(30.0, 0.0, 1.0, 0.5);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_y_equals_mu_smooth) {
  using generalized_normal_lpdf_test::expect_all_rotations;

  // expect_ad checks up to third derivatives, and d3/dy3 log p only exists
  // at y == mu for beta > 3, so the tie is only tested there.
  expect_all_rotations(0.5, 0.5, 1.3, 4.0);
  expect_all_rotations(-2.0, -2.0, 0.4, 4.5);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_y_equals_mu_normal) {
  using generalized_normal_lpdf_test::f_y_mu_alpha;

  // With beta held fixed at 2 the distribution is normal(mu,
  // alpha / sqrt(2)), so all y, mu, and alpha derivatives exist.
  stan::test::expect_ad(f_y_mu_alpha(2.0), 0.5, 0.5, 1.3);
  stan::test::expect_ad(f_y_mu_alpha(2.0), 0.0, 0.0, 1.0);

  // Same requirement inside a vectorized call, where only y[0] == mu[0].
  auto f_vec = [](const auto& y, const auto& mu, const auto& alpha) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, 2.0);
  };
  Eigen::VectorXd y(2), mu(2), alpha(2);
  y << 0.5, 0.3;
  mu << 0.5, -0.2;
  alpha << 1.0, 1.4;
  stan::test::expect_ad(f_vec, y, mu, alpha);

  // The same boundary handling applies elementwise when beta is a vector.
  Eigen::VectorXd beta(2);
  beta << 2.0, 3.5;
  auto f_vec_beta = [&beta](const auto& y, const auto& mu, const auto& alpha) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, beta);
  };
  stan::test::expect_ad(f_vec_beta, y, mu, alpha);
}

TEST_F(AgradRev,
       mathMixScalFun_generalized_normal_lpdf_y_equals_mu_variable_beta) {
  using stan::math::generalized_normal_lpdf;

  // Pure beta derivatives are smooth at the tie. In particular, the kernel
  // contribution r^beta * log(r) is zero by continuity.
  auto f_beta = [](const auto& beta) {
    return generalized_normal_lpdf(0.5, 0.5, 1.3, beta);
  };
  stan::test::expect_ad(f_beta, 2.0);

  // The joint value, gradient, and Hessian also exist at beta == 2. Do not use
  // expect_ad on all four arguments here: it additionally checks the mixed
  // third derivative d3/(d beta d y2), which diverges at y == mu.
  auto f_joint = [](const auto& x) {
    return generalized_normal_lpdf(x(0), x(1), x(2), x(3));
  };
  Eigen::VectorXd x(4);
  x << 0.5, 0.5, 1.3, 2.0;
  double fx;
  Eigen::VectorXd grad;
  Eigen::MatrixXd hessian;
  stan::math::hessian(f_joint, x, fx, grad, hessian);

  const double alpha = x(2);
  EXPECT_NEAR(-stan::math::LOG_TWO - std::log(alpha) - stan::math::lgamma(1.5),
              fx, 1e-12);
  EXPECT_NEAR(0.0, grad(0), 1e-12);
  EXPECT_NEAR(0.0, grad(1), 1e-12);
  EXPECT_NEAR(-1.0 / alpha, grad(2), 1e-12);
  EXPECT_NEAR(stan::math::digamma(1.5) / 4.0, grad(3), 1e-12);

  const double normal_curvature = -2.0 / (alpha * alpha);
  EXPECT_NEAR(normal_curvature, hessian(0, 0), 1e-12);
  EXPECT_NEAR(-normal_curvature, hessian(0, 1), 1e-12);
  EXPECT_NEAR(-normal_curvature, hessian(1, 0), 1e-12);
  EXPECT_NEAR(normal_curvature, hessian(1, 1), 1e-12);
  EXPECT_NEAR(0.0, hessian(0, 3), 1e-12);
  EXPECT_NEAR(0.0, hessian(1, 3), 1e-12);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_y_equals_mu_cusp) {
  using stan::math::generalized_normal_lpdf;
  using stan::math::var;

  // At the tie, the first derivative is undefined for beta <= 1 and the
  // second derivative is not finite for beta < 2. Check Stan's zero-gradient
  // convention and the remaining parameter derivatives directly.
  for (double beta : {0.5, 1.0, 1.5}) {
    var y = 0.75;
    var mu = 0.75;
    var alpha = 1.4;
    var beta_v = beta;
    var lp = generalized_normal_lpdf(y, mu, alpha, beta_v);

    double expect_lp = -stan::math::LOG_TWO - std::log(1.4)
                       - stan::math::lgamma(1.0 + 1.0 / beta);
    EXPECT_FLOAT_EQ(expect_lp, lp.val());

    lp.grad();
    EXPECT_FLOAT_EQ(0.0, y.adj());
    EXPECT_FLOAT_EQ(0.0, mu.adj());
    EXPECT_FLOAT_EQ(-1.0 / 1.4, alpha.adj());
    EXPECT_FLOAT_EQ(stan::math::digamma(1.0 + 1.0 / beta) / (beta * beta),
                    beta_v.adj());
    stan::math::recover_memory();
  }
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_vectorized) {
  auto f_vec = [](const auto& y, const auto& mu, const auto& alpha) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, 3.0);
  };
  auto f_vec_beta = [](const auto& y, const auto& mu, const auto& beta) {
    return stan::math::generalized_normal_lpdf(y, mu, 1.1, beta);
  };

  Eigen::VectorXd y(3), mu(3), alpha(3), beta(3);
  y << -1.0, 0.6, 2.0;
  mu << 0.0, 0.5, -0.4;
  alpha << 1.0, 2.5, 0.7;
  beta << 3.5, 4.0, 5.5;

  stan::test::expect_ad(f_vec, y, mu, alpha);
  stan::test::expect_ad(f_vec, y, mu[0], alpha);
  stan::test::expect_ad(f_vec, y, mu, alpha[0]);
  stan::test::expect_ad(f_vec, y[0], mu, alpha);
  stan::test::expect_ad(f_vec, y[1], mu[0], alpha);

  stan::test::expect_ad(f_vec_beta, y, mu, beta);
  stan::test::expect_ad(f_vec_beta, y, mu, beta[0]);
  stan::test::expect_ad(f_vec_beta, y[1], mu[0], beta);

  std::vector<double> y_std{-1.0, 0.6, 2.0};
  std::vector<double> mu_std{0.0, 0.5, -0.4};
  std::vector<double> alpha_std{1.0, 2.5, 0.7};
  stan::test::expect_ad(f_vec, y_std, mu_std, alpha_std);
  stan::test::expect_ad(f_vec, y_std, mu_std[0], alpha_std);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_vectorized_tie) {
  // y[1] == mu[1] exercises the y == mu branch inside a vectorized call;
  // beta > 3 keeps the third derivative at the tie well defined.
  auto f_vec = [](const auto& y, const auto& mu, const auto& alpha) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, 4.0);
  };
  auto f_vec_beta = [](const auto& y, const auto& mu, const auto& beta) {
    return stan::math::generalized_normal_lpdf(y, mu, 1.1, beta);
  };

  Eigen::VectorXd y(3), mu(3), alpha(3), beta(3);
  y << -1.0, 0.5, 2.0;
  mu << 0.0, 0.5, -0.4;
  alpha << 1.0, 2.5, 0.7;
  beta << 3.5, 4.0, 5.5;

  stan::test::expect_ad(f_vec, y, mu, alpha);
  stan::test::expect_ad(f_vec, y, mu, alpha[0]);
  stan::test::expect_ad(f_vec, y[1], mu[1], alpha);
  stan::test::expect_ad(f_vec_beta, y, mu, beta);
  stan::test::expect_ad(f_vec_beta, y[1], mu[1], beta);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_propto) {
  using stan::math::generalized_normal_lpdf;
  using stan::math::var;

  // expect_ad cannot be used here: with all-double arguments propto=true
  // drops every term and returns 0, so no value comparison is possible.
  for (double y_val : {1.3, 0.4}) {
    Eigen::Matrix<var, -1, 1> x(4);
    x << y_val, 0.4, 1.7, 2.5;
    var lp_propto = generalized_normal_lpdf<true>(x(0), x(1), x(2), x(3));
    lp_propto.grad();
    Eigen::VectorXd grad_propto(4);
    for (int i = 0; i < 4; ++i)
      grad_propto(i) = x(i).adj();
    double val_propto = lp_propto.val();
    stan::math::recover_memory();

    Eigen::Matrix<var, -1, 1> x2(4);
    x2 << y_val, 0.4, 1.7, 2.5;
    var lp_full = generalized_normal_lpdf<false>(x2(0), x2(1), x2(2), x2(3));
    lp_full.grad();

    // Only the -log(2) normalizing constant is dropped when propto=true.
    EXPECT_FLOAT_EQ(lp_full.val() + stan::math::LOG_TWO, val_propto);
    for (int i = 0; i < 4; ++i)
      EXPECT_FLOAT_EQ(x2(i).adj(), grad_propto(i));
    stan::math::recover_memory();
  }
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_errors) {
  using generalized_normal_lpdf_test::f_mu_alpha_beta;
  using generalized_normal_lpdf_test::f_y_alpha_beta;
  using generalized_normal_lpdf_test::f_y_mu_alpha;
  using generalized_normal_lpdf_test::f_y_mu_beta;

  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();

  // alpha must be positive, beta must be positive, mu must be finite,
  // y must not be NaN -- every instantiation must throw where double does.
  stan::test::expect_ad(f_y_mu_alpha(2.0), 1.0, 0.0, 0.0);
  stan::test::expect_ad(f_y_mu_alpha(2.0), 1.0, 0.0, -1.5);
  stan::test::expect_ad(f_y_mu_beta(1.0), 1.0, 0.0, 0.0);
  stan::test::expect_ad(f_y_mu_beta(1.0), 1.0, 0.0, -2.0);
  stan::test::expect_ad(f_y_alpha_beta(0.0), nan, 1.0, 2.0);
  stan::test::expect_ad(f_mu_alpha_beta(1.0), inf, 1.0, 2.0);
  stan::test::expect_ad(f_mu_alpha_beta(1.0), -inf, 1.0, 2.0);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_infinite_params) {
  using stan::math::generalized_normal_lpdf;
  double inf = std::numeric_limits<double>::infinity();

  // Infinite scale and shape parameters are outside the supported domain.
  EXPECT_THROW(generalized_normal_lpdf(1.0, 0.0, inf, 2.0), std::domain_error);
  EXPECT_THROW(generalized_normal_lpdf(1.0, 0.0, 1.0, inf), std::domain_error);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_subnormal_alpha) {
  using stan::math::generalized_normal_lpdf;

  // inv(alpha) overflows for subnormal alpha, so abs(diff) * inv_alpha is
  // 0 * inf = NaN at y == mu and the |y - mu| == 0 test never fires.
  for (double a : {1e-309, 1e-310, std::numeric_limits<double>::denorm_min()}) {
    double expected
        = -stan::math::LOG_TWO - std::log(a) - stan::math::lgamma(1.5);
    EXPECT_FLOAT_EQ(expected, generalized_normal_lpdf(0.0, 0.0, a, 2.0));
  }
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_matvar) {
  auto f = [](const auto& y, const auto& mu, const auto& alpha) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, 3.0);
  };
  auto f_beta = [](const auto& y, const auto& mu, const auto& beta) {
    return stan::math::generalized_normal_lpdf(y, mu, 1.1, beta);
  };

  Eigen::VectorXd y(3), mu(3), alpha(3), beta(3);
  y << -1.0, 0.6, 2.0;
  mu << 0.0, 0.5, -0.4;
  alpha << 1.0, 2.5, 0.7;
  beta << 3.5, 4.0, 5.5;

  stan::test::expect_ad_matvar(f, y, mu, alpha);
  stan::test::expect_ad_matvar(f, y, mu, alpha[0]);
  stan::test::expect_ad_matvar(f_beta, y, mu, beta);
}

TEST_F(AgradRev, mathMixScalFun_generalized_normal_lpdf_size_errors) {
  using stan::math::generalized_normal_lpdf;

  Eigen::VectorXd y(3), mu(2);
  y << 0.1, 0.2, 0.3;
  mu << 0.0, 0.0;
  EXPECT_THROW(generalized_normal_lpdf(y, mu, 1.0, 2.0), std::invalid_argument);

  Eigen::VectorXd empty(0);
  EXPECT_FLOAT_EQ(0.0, generalized_normal_lpdf(empty, 0.0, 1.0, 2.0));
}
