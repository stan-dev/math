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

  beta(1) = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(categorical_logit_rng(beta, rng));

  beta(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(categorical_logit_rng(beta, rng), std::domain_error);

  beta << -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity();
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

TEST(ProbDistributionsCategoricalLogit, multiplePosInfinityIsUniform) {
  using Eigen::VectorXd;
  using stan::math::categorical_logit_rng;
  boost::random::mt19937 rng;

  VectorXd beta(4);
  beta << -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), 2.0,
      std::numeric_limits<double>::infinity();

  int N = 10000;
  int count_2 = 0;
  int count_4 = 0;
  for (int i = 0; i < N; i++) {
    int result = categorical_logit_rng(beta, rng);
    if (result == 2) count_2++;
    if (result == 4) count_4++;
  }

  // every draw lands on one of the +inf entries, split roughly evenly
  EXPECT_EQ(N, count_2 + count_4);
  EXPECT_NEAR(count_2, N / 2, N * 0.05);
}
