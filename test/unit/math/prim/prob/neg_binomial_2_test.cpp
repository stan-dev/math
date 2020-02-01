#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/NegativeBinomial2LogTestRig.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>
#include <string>

class NegativeBinomial2TestRig : public VectorIntRNGTestRig {
 public:
  NegativeBinomial2TestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {0.1, 1.7, 3.99},
                            {1, 2, 3}, {-2.1, -0.5, 0.0}, {-3, -1, 0},
                            {0.1, 1.1, 4.99}, {1, 2, 3}, {-3.0, -2.0, 0.0},
                            {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& phi, const T3&,
                        T_rng& rng) const {
    return stan::math::neg_binomial_2_rng(mu, phi, rng);
  }

  template <typename T1>
  double pmf(int y, T1 mu, double phi, double) const {
    return std::exp(stan::math::neg_binomial_2_lpmf(y, mu, phi));
  }
};

TEST(ProbDistributionsNegativeBinomial2, errorCheck) {
  check_dist_throws_all_types(NegativeBinomial2TestRig());
}

TEST(ProbDistributionsNegativeBinomial2, distributionCheck) {
  check_counts_real_real(NegativeBinomial2TestRig());
}

TEST(ProbDistributionsNegBinomial2, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::neg_binomial_2_rng(6, 2, rng));
  EXPECT_NO_THROW(stan::math::neg_binomial_2_rng(0.5, 1, rng));
  EXPECT_NO_THROW(stan::math::neg_binomial_2_rng(1e8, 1, rng));

  EXPECT_THROW(stan::math::neg_binomial_2_rng(0, -2, rng), std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_rng(6, -2, rng), std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_rng(-6, -0.1, rng),
               std::domain_error);
  EXPECT_THROW(
      stan::math::neg_binomial_2_rng(stan::math::positive_infinity(), 2, rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::neg_binomial_2_rng(stan::math::positive_infinity(), 6, rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::neg_binomial_2_rng(2, stan::math::positive_infinity(), rng),
      std::domain_error);

  std::string error_msg;

  error_msg
      = "neg_binomial_2_rng: Location parameter "
        "divided by the precision parameter is "
        "inf, but must be finite!";
  try {
    stan::math::neg_binomial_2_rng(1e300, 1e-300, rng);
    FAIL() << "neg_binomial_2_rng should have thrown" << std::endl;
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(error_msg) == std::string::npos)
      FAIL() << "Error message is different than expected" << std::endl
             << "EXPECTED: " << error_msg << std::endl
             << "FOUND: " << e.what() << std::endl;
    SUCCEED();
  }

  error_msg
      = "neg_binomial_2_rng: Random number that "
        "came from gamma distribution is";
  try {
    stan::math::neg_binomial_2_rng(1e10, 1e20, rng);
    FAIL() << "neg_binomial_2_rng should have thrown" << std::endl;
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(error_msg) == std::string::npos)
      FAIL() << "Error message is different than expected" << std::endl
             << "EXPECTED: " << error_msg << std::endl
             << "FOUND: " << e.what() << std::endl;
    SUCCEED();
  }
}

TEST(ProbDistributionsNegBinomial2, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::negative_binomial_distribution<> dist(1.1, 1.1 / (1.1 + 2.4));
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
    int a = stan::math::neg_binomial_2_rng(2.4, 1.1, rng);
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

TEST(ProbDistributionsNegBinomial2, chiSquareGoodnessFitTest2) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::negative_binomial_distribution<> dist(0.6, 0.6 / (0.6 + 2.4));
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
    int a = stan::math::neg_binomial_2_rng(2.4, 0.6, rng);
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

TEST(ProbDistributionsNegBinomial2, chiSquareGoodnessFitTest3) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::negative_binomial_distribution<> dist(30, 30 / (30 + 60.4));
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
    int a = stan::math::neg_binomial_2_rng(60.4, 30, rng);
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

TEST(ProbDistributionsNegBinomial, chiSquareGoodnessFitTest4) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::negative_binomial_distribution<> dist(80, 80 / (80 + 30.4));
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
    int a = stan::math::neg_binomial_2_rng(30.4, 80, rng);
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

TEST(ProbDistributionsNegBinomial, extreme_values) {
  int N = 100;
  double mu = 8;
  double phi = 1e12;
  for (int n = 0; n < 10; ++n) {
    phi *= 10;
    double logp = stan::math::neg_binomial_2_log<false>(N, mu, phi);
    EXPECT_LT(logp, 0);
  }
}

TEST(ProbDistributionsNegBinomial2, vectorAroundCutoff) {
  int y = 10;
  double mu = 9.36;
  std::vector<double> phi;
  phi.push_back(1);
  phi.push_back(1e15);
  double vector_value = stan::math::neg_binomial_2_lpmf(y, mu, phi);
  double scalar_value = stan::math::neg_binomial_2_lpmf(y, mu, phi[0])
                        + stan::math::neg_binomial_2_lpmf(y, mu, phi[1]);

  EXPECT_FLOAT_EQ(vector_value, scalar_value);
}

TEST(ProbDistributionsNegativeBinomial2Log, distributionCheck) {
  check_counts_real_real(NegativeBinomial2LogTestRig());
}
