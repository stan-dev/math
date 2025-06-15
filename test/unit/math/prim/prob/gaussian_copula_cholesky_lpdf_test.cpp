#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDistributionsGaussianCopulaCholesky, EqualsMultiNormalCholesky) {
  using stan::math::gaussian_copula_cholesky_lpdf;
  using stan::math::multi_normal_cholesky_lpdf;
  using stan::math::diag_pre_multiply;

  Eigen::VectorXd y(2);
  y << 1, 3;

  Eigen::VectorXd mu(2);
  mu << 0.1, 0;

  Eigen::VectorXd sd(2);
  sd << 2, 1;

  Eigen::MatrixXd chol(2, 2);
  chol << 1, 0, 0.5, 0.8660254037844386;

  auto norm_lcdf = [](auto&& y, auto&& mu, auto&& sigma) {
    return stan::math::normal_lcdf(y, mu, sigma);
  };

  auto std_norm_lcdf = [](auto&& y) {
    return stan::math::std_normal_lcdf(y);
  };

  // y[0] ~ normal(0.1, 2)
  // y[1] ~ normal(0, 1)
  auto lcdf_functors = std::make_tuple(
    std::make_tuple(norm_lcdf, mu[0], sd[0]),
    std::make_tuple(std_norm_lcdf)
  );

  double log_prob =
    stan::math::normal_lpdf(y[0], 0.1, 2) +
    stan::math::std_normal_lpdf(y[1]) +
    gaussian_copula_cholesky_lpdf(y, lcdf_functors, chol);


  double expected_log_prob = multi_normal_cholesky_lpdf(
    y, mu, diag_pre_multiply(sd, chol));

  EXPECT_FLOAT_EQ(log_prob, expected_log_prob);
}

TEST(ProbDistributionsGaussianCopulaCholesky, NonNormalMarginals) {
  Eigen::VectorXd y(2);
  y << 10, 6;

  Eigen::MatrixXd chol(2, 2);
  chol << 1, 0, 0.2, 0.9797958971;

  auto gamma_lcdf = [](auto&& y, auto&& shape, auto&& scale) {
    return stan::math::gamma_lcdf(y, shape, scale);
  };

  auto exp_lcdf = [](auto&& y, auto&& scale) {
    return stan::math::exponential_lcdf(y, scale);
  };

  // y[0] ~ gamma(2, 1)
  // y[1] ~ exponential(2)
  auto lcdf_functors = std::make_tuple(
    std::make_tuple(gamma_lcdf, 2.0, 1.0),
    std::make_tuple(exp_lcdf, 2.0)
  );

  double log_prob =
    stan::math::gamma_lpdf(y[0], 2.0, 1.0) +
    stan::math::exponential_lpdf(y[1], 2.0) +
    stan::math::gaussian_copula_cholesky_lpdf(y, lcdf_functors, chol);

  EXPECT_FLOAT_EQ(log_prob, -16.61005941);
}

TEST(ProbDistributionsGaussianCopulaCholesky, Errors) {
  Eigen::VectorXd y(2);
  y << 10, 6;

  Eigen::VectorXd small_y(1);
  small_y << 10;

  Eigen::MatrixXd chol(2, 2);
  chol << 1, 0, 0.2, 0.9797958971;

  auto gamma_lcdf = [](auto&& y, auto&& shape, auto&& scale) {
    return stan::math::gamma_lcdf(y, shape, scale);
  };
  auto exp_lcdf = [](auto&& y, auto&& scale) {
    return stan::math::exponential_lcdf(y, scale);
  };

  auto lcdf_functors = std::make_tuple(
    std::make_tuple(gamma_lcdf, 2.0, 1.0),
    std::make_tuple(exp_lcdf, 2.0)
  );

  auto small_lcdf_functors = std::make_tuple(
    std::make_tuple(gamma_lcdf, 2.0, 1.0)
  );

  auto invalid_lcdf_functors = std::make_tuple(
    std::make_tuple(gamma_lcdf, 2.0, 1.0),
    std::make_tuple([](auto&& y){ return y * 10; })
  );

  EXPECT_THROW(
    stan::math::gaussian_copula_cholesky_lpdf(y, small_lcdf_functors, chol),
    std::invalid_argument
  );

  EXPECT_THROW(
    stan::math::gaussian_copula_cholesky_lpdf(small_y, lcdf_functors, chol),
    std::invalid_argument
  );

  EXPECT_THROW(
    stan::math::gaussian_copula_cholesky_lpdf(y, invalid_lcdf_functors, chol),
    std::domain_error
  );
}
