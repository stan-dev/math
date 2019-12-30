#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/scal.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>
#include <string>

class NegativeBinomialTestRig : public VectorIntRNGTestRig {
 public:
  NegativeBinomialTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {0.1, 1.7, 3.99},
                            {1, 2, 3}, {-2.1, -0.5, 0.0}, {-3, -1, 0},
                            {0.1, 1.1, 4.99}, {1, 2, 3}, {-3.0, -2.0, 0.0},
                            {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& alpha, const T2& beta, const T3&,
                        T_rng& rng) const {
    return stan::math::neg_binomial_rng(alpha, beta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 alpha, double beta, double) const {
    return std::exp(stan::math::neg_binomial_lpmf(y, alpha, beta));
  }
};

TEST(ProbDistributionsNegativeBinomial, errorCheck) {
  check_dist_throws_all_types(NegativeBinomialTestRig());
}

TEST(ProbDistributionsNegativeBinomial, distributionCheck) {
  check_counts_real_real(NegativeBinomialTestRig());
}

TEST(ProbDistributionsNegBinomial, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::neg_binomial_rng(6, 2, rng));
  EXPECT_NO_THROW(stan::math::neg_binomial_rng(0.5, 1, rng));
  EXPECT_NO_THROW(stan::math::neg_binomial_rng(1e9, 1, rng));

  EXPECT_THROW(stan::math::neg_binomial_rng(0, -2, rng), std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_rng(6, -2, rng), std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_rng(-6, -0.1, rng), std::domain_error);
  EXPECT_THROW(
      stan::math::neg_binomial_rng(stan::math::positive_infinity(), 2, rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::neg_binomial_rng(stan::math::positive_infinity(), 6, rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::neg_binomial_rng(2, stan::math::positive_infinity(), rng),
      std::domain_error);

  std::string error_msg;
  error_msg
      = "neg_binomial_rng: Random number that "
        "came from gamma distribution is";
  try {
    stan::math::neg_binomial_rng(1e10, 1, rng);
    FAIL() << "neg_binomial_rng should have thrown" << std::endl;
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(error_msg) == std::string::npos)
      FAIL() << "Error message is different than expected" << std::endl
             << "EXPECTED: " << error_msg << std::endl
             << "FOUND: " << e.what() << std::endl;
    SUCCEED();
  }
}

void expected_bin_sizes(double* expect, const int K, const int N,
                        const double alpha, const double beta) {
  double p = 0;
  for (int i = 0; i < K; i++) {
    expect[i] = N * std::exp(stan::math::neg_binomial_log(i, alpha, beta));
    p += std::exp(stan::math::neg_binomial_log(i, alpha, beta));
  }
  expect[K - 1] = N * (1.0 - p);
}

TEST(ProbDistributionsNegBinomial, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  double p = 0.6;
  double alpha = 5;
  double beta = p / (1 - p);
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, (1 - p)));
  boost::math::chi_squared mydist(K - 1);

  int loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = i - 1;

  int count = 0;
  double bin[K];
  double expect[K];

  for (int i = 0; i < K; i++)
    bin[i] = 0;
  expected_bin_sizes(expect, K, N, alpha, beta);

  while (count < N) {
    int a = stan::math::neg_binomial_rng(alpha, beta, rng);
    int i = 0;
    while (i < K - 1 && a > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_LT(chi, boost::math::quantile(boost::math::complement(mydist, 1e-6)));
}

TEST(ProbDistributionsNegBinomial, chiSquareGoodnessFitTest2) {
  boost::random::mt19937 rng;
  double p = 0.8;
  double alpha = 2.4;
  double beta = p / (1 - p);
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, (1 - p)));
  boost::math::chi_squared mydist(K - 1);

  int loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = i - 1;

  int count = 0;
  double bin[K];
  double expect[K];

  for (int i = 0; i < K; i++)
    bin[i] = 0;
  expected_bin_sizes(expect, K, N, alpha, beta);

  while (count < N) {
    int a = stan::math::neg_binomial_rng(alpha, beta, rng);
    int i = 0;
    while (i < K - 1 && a > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_LT(chi, boost::math::quantile(boost::math::complement(mydist, 1e-6)));
}

TEST(ProbDistributionsNegBinomial, chiSquareGoodnessFitTest3) {
  boost::random::mt19937 rng;
  double p = 0.2;
  double alpha = 0.4;
  double beta = p / (1 - p);
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, (1 - p)));
  boost::math::chi_squared mydist(K - 1);

  int loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = i - 1;

  int count = 0;
  double bin[K];
  double expect[K];

  for (int i = 0; i < K; i++)
    bin[i] = 0;
  expected_bin_sizes(expect, K, N, alpha, beta);

  while (count < N) {
    int a = stan::math::neg_binomial_rng(alpha, beta, rng);
    int i = 0;
    while (i < K - 1 && a > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_LT(chi, boost::math::quantile(boost::math::complement(mydist, 1e-6)));
}
