#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

class PoissonTestRig : public VectorIntRNGTestRig {
 public:
  PoissonTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {0.1, 1.1, 4.99},
                            {1, 2, 3}, {-3.0, -2.0, 0.0}, {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& lambda, const T2&, const T3&,
                        T_rng& rng) const {
    return stan::math::poisson_rng(lambda, rng);
  }

  template <typename T1>
  double pmf(int y, T1 lambda, double, double) const {
    return std::exp(stan::math::poisson_lpmf(y, lambda));
  }
};

TEST(ProbDistributionsPoisson, errorCheck) {
  check_dist_throws_all_types(PoissonTestRig());
}

TEST(ProbDistributionsPoisson, distributionCheck) {
  check_counts_real(PoissonTestRig());
}

TEST(ProbDistributionsPoisson, error_check) {
  using std::log;

  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::poisson_rng(6, rng));

  EXPECT_THROW(stan::math::poisson_rng(-6, rng), std::domain_error);

  EXPECT_NO_THROW(stan::math::poisson_rng(1e9, rng));

  EXPECT_THROW(stan::math::poisson_rng(pow(2.0, 31), rng), std::domain_error);

  EXPECT_NO_THROW(stan::math::poisson_log_rng(6, rng));

  EXPECT_NO_THROW(stan::math::poisson_log_rng(-6, rng));

  EXPECT_NO_THROW(stan::math::poisson_log_rng(log(1e9), rng));

  EXPECT_THROW(stan::math::poisson_log_rng(log(pow(2.0, 31)), rng),
               std::domain_error);
}

TEST(ProbDistributionsPoisson, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::poisson_distribution<> dist(5);
  boost::math::chi_squared mydist(K - 1);

  int loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = i - 1;

  int count = 0;
  double bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N * pdf(dist, i);
  }
  expect[K - 1] = N * (1 - cdf(dist, K - 1));

  while (count < N) {
    int a = stan::math::poisson_rng(5, rng);
    int i = 0;
    while (i < K - 1 && a > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsPoisson, chiSquareGoodnessFitTest2) {
  using std::log;

  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::poisson_distribution<> dist(5);
  boost::math::chi_squared mydist(K - 1);

  int loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = i - 1;

  int count = 0;
  double bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N * pdf(dist, i);
  }
  expect[K - 1] = N * (1 - cdf(dist, K - 1));

  while (count < N) {
    int a = stan::math::poisson_log_rng(log(5), rng);
    int i = 0;
    while (i < K - 1 && a > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbPoisson, log_matches_lpmf) {
  int y = 3;
  double lambda = 2.3;

  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf(y, lambda)),
                  (stan::math::poisson_log(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<true>(y, lambda)),
                  (stan::math::poisson_log<true>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<false>(y, lambda)),
                  (stan::math::poisson_log<false>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<true, int, double>(y, lambda)),
                  (stan::math::poisson_log<true, int, double>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<false, int, double>(y, lambda)),
                  (stan::math::poisson_log<false, int, double>(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lpmf<int, double>(y, lambda)),
                  (stan::math::poisson_log<int, double>(y, lambda)));
}
