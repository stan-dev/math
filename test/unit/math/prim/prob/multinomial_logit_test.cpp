#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbDistributionsMultinomialLogit, RNGZero) {
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> beta(3);
  beta << 1.3, 0.1, -2.6;
  // bug in 4.8.1: RNG does not allow a zero total count
  EXPECT_NO_THROW(stan::math::multinomial_logit_rng(beta, 0, rng));
  // when the total count is zero, the sample should be a zero array
  std::vector<int> sample = stan::math::multinomial_logit_rng(beta, 0, rng);
  for (int k : sample) {
    EXPECT_EQ(0, k);
  }
}

TEST(ProbDistributionsMultinomialLogit, RNGSize) {
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> beta(5);
  beta << log(0.3), log(0.1), log(0.2), log(0.2), log(0.2);
  std::vector<int> sample = stan::math::multinomial_logit_rng(beta, 10, rng);
  EXPECT_EQ(5U, sample.size());
}

TEST(ProbDistributionsMultinomialLogit, RNGErrorCheck) {
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> beta(3);

  beta << 1.3, 0.1, -2.6;
  EXPECT_NO_THROW(stan::math::multinomial_logit_rng(beta, 10, rng));

  // +inf allowed
  beta(1) = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(stan::math::multinomial_logit_rng(beta, 10, rng));

  // NaN throws
  beta(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::multinomial_logit_rng(beta, 10, rng),
               std::domain_error);

  // all -inf throws
  beta << -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multinomial_logit_rng(beta, 10, rng),
               std::domain_error);
}

TEST(ProbDistributionsMultinomialLogit, RNGMultiplePosInfinityIsUniform) {
  boost::random::mt19937 rng;
  // multiple +inf case: uniform over the +inf entries
  Matrix<double, Dynamic, 1> beta(4);
  beta << -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), 2.0,
      std::numeric_limits<double>::infinity();

  int N = 20;
  int trials = 1000;
  int total_2 = 0;

  for (int i = 0; i < trials; i++) {
    std::vector<int> sample = stan::math::multinomial_logit_rng(beta, N, rng);
    // -inf and finite entries get no counts
    EXPECT_EQ(0, sample[0]);
    EXPECT_EQ(0, sample[2]);
    EXPECT_EQ(N, sample[1] + sample[3]);
    total_2 += sample[1];
  }

  // roughly even split of all counts between the two +inf entries
  double expected = trials * N / 2.0;
  EXPECT_NEAR(total_2, expected, expected * 0.05);
}

TEST(ProbDistributionsMultinomialLogit, MultinomialLogit) {
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);
  EXPECT_FLOAT_EQ(-2.002481, stan::math::multinomial_logit_lpmf(ns, beta));
}
TEST(ProbDistributionsMultinomialLogit, Propto) {
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);
  EXPECT_FLOAT_EQ(0.0, stan::math::multinomial_logit_lpmf<true>(ns, beta));
}

using stan::math::multinomial_logit_lpmf;

TEST(ProbDistributionsMultinomialLogit, error) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);

  EXPECT_NO_THROW(multinomial_logit_lpmf(ns, beta));

  ns[1] = 0;
  EXPECT_NO_THROW(multinomial_logit_lpmf(ns, beta));
  ns[1] = -1;
  EXPECT_THROW(multinomial_logit_lpmf(ns, beta), std::domain_error);
  ns[1] = 1;

  beta(0) = nan;
  EXPECT_THROW(multinomial_logit_lpmf(ns, beta), std::domain_error);
  // +/-inf allowed for data
  beta(0) = inf;
  EXPECT_NO_THROW(multinomial_logit_lpmf(ns, beta));
  beta(0) = -inf;
  EXPECT_NO_THROW(multinomial_logit_lpmf(ns, beta));
  // all -inf throws
  beta << -inf, -inf, -inf;
  EXPECT_THROW(multinomial_logit_lpmf(ns, beta), std::domain_error);

  beta(0) = 0.2;
  beta(1) = 0.3;
  beta(2) = 0.5;

  ns.resize(2);
  EXPECT_THROW(multinomial_logit_lpmf(ns, beta), std::invalid_argument);
}

TEST(ProbDistributionsMultinomialLogit, zeros) {
  double result;
  std::vector<int> ns;
  ns.push_back(0);
  ns.push_back(1);
  ns.push_back(2);
  Matrix<double, Dynamic, 1> beta(3, 1);
  beta << log(0.2), log(0.3), log(0.5);

  result = multinomial_logit_lpmf(ns, beta);
  EXPECT_FALSE(std::isnan(result));

  std::vector<int> ns2;
  ns2.push_back(0);
  ns2.push_back(0);
  ns2.push_back(0);

  double result2 = multinomial_logit_lpmf(ns2, beta);
  EXPECT_FLOAT_EQ(0.0, result2);
}

TEST(ProbDistributionsMultinomialLogit, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int M = 10;
  int trials = 1000;
  int N = M * trials;

  int K = 3;
  Matrix<double, Dynamic, 1> beta(K);
  beta << -1, 1, -10;
  Eigen::VectorXd theta = stan::math::softmax(beta);
  boost::math::chi_squared mydist(K - 1);

  double expect[K];
  for (int i = 0; i < K; ++i)
    expect[i] = N * theta(i);

  int bin[K];
  for (int i = 0; i < K; ++i)
    bin[i] = 0;

  for (int count = 0; count < M; ++count) {
    std::vector<int> a = stan::math::multinomial_logit_rng(beta, trials, rng);
    for (int i = 0; i < K; ++i)
      bin[i] += a[i];
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j])) / expect[j];

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsMultinomialLogit, infinityValues) {
  double inf = std::numeric_limits<double>::infinity();

  // -inf entries have zero probability, others as usual
  Matrix<double, Dynamic, 1> beta(3);
  beta << -inf, 1.0, 2.0;
  Matrix<double, Dynamic, 1> finite(2);
  finite << 1.0, 2.0;
  std::vector<int> ns{0, 1, 2};
  std::vector<int> finite_ns{1, 2};
  EXPECT_FLOAT_EQ(multinomial_logit_lpmf(finite_ns, finite),
                  multinomial_logit_lpmf(ns, beta));
  // nonzero count on a -inf entry is impossible
  ns = {1, 1, 1};
  EXPECT_EQ(-inf, multinomial_logit_lpmf(ns, beta));

  // +inf entries split the probability evenly
  beta.resize(4);
  beta << -inf, inf, 2.0, inf;
  ns = {0, 1, 0, 2};
  EXPECT_FLOAT_EQ(std::log(3.0) - 3 * std::log(2.0),
                  multinomial_logit_lpmf(ns, beta));
  // nonzero count outside the +inf entries is impossible
  ns = {0, 1, 1, 1};
  EXPECT_EQ(-inf, multinomial_logit_lpmf(ns, beta));
}
