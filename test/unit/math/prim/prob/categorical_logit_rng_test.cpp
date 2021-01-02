#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>

TEST(ProbDistributionsCategoricalLogit, error_check) {
  using Eigen::VectorXd;
  using stan::math::categorical_logit_rng;
  boost::random::mt19937 rng;

  VectorXd beta(3);

  beta << 1.0, 10.0, -10.0;
  EXPECT_NO_THROW(categorical_logit_rng(beta, rng));

  beta << -1e3, 1.1e3, 1e5;
  EXPECT_NO_THROW(categorical_logit_rng(beta, rng));

  beta(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(categorical_logit_rng(beta, rng), std::domain_error);

  beta(1) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(categorical_logit_rng(beta, rng), std::domain_error);
}

TEST(ProbDistributionsCategoricalLogit, chiSquareGoodnessFitTest) {
  using Eigen::VectorXd;
  using stan::math::softmax;
  boost::random::mt19937 rng;
  int N = 10000;
  int K = 3;
  VectorXd beta(K);

  beta << -0.5, 0.1, 0.3;

  VectorXd theta = softmax(beta);
  boost::math::chi_squared mydist(K - 1);

  int bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N * theta(i);
  }

  for (int i = 0; i < N; i++) {
    int a = stan::math::categorical_logit_rng(beta, rng);
    bin[a - 1]++;
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
