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

TEST(ProbDistributionsCategoricalLogit, categorical_logit_vecRNG) {
  using stan::math::categorical_logit_rng;
  boost::random::mt19937 rng;

  Eigen::VectorXd vec1(4);
  vec1 << Eigen::VectorXd::Random(4);
  Eigen::VectorXd vec2(4);
  vec2 << Eigen::VectorXd::Random(4);
  Eigen::VectorXd vec3(4);
  vec3 << Eigen::VectorXd::Random(4);

  std::vector<Eigen::VectorXd> vecs{vec1, vec2, vec3};

  std::vector<int> rng_stvec = categorical_logit_rng(vecs, rng);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(rng_stvec[i], categorical_logit_rng(vecs[i], rng));
  }
}

TEST(ProbDistributionsCategoricalLogit, categorical_logit_vecRNG_throw) {
  using stan::math::categorical_logit_rng;
  using stan::math::INFTY;
  using stan::math::softmax;
  boost::random::mt19937 rng;

  Eigen::VectorXd vec1(4);
  vec1 << Eigen::VectorXd::Random(4);
  Eigen::VectorXd vec2(4);
  vec2 << 0.5, -1.2, INFTY, 1.5;
  Eigen::VectorXd vec3(4);
  vec3 << Eigen::VectorXd::Random(4);

  std::vector<Eigen::VectorXd> vecs{vec1, vec2, vec3};

  EXPECT_THROW(categorical_logit_rng(vecs, rng), std::domain_error);
}
