#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/NegativeBinomial2LogTestRig.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>

TEST(ProbDistributionsNegativeBinomial2Log, errorCheck) {
  check_dist_throws_all_types(NegativeBinomial2LogTestRig());
}

TEST(ProbDistributionsNegBinomial2Log, error_check) {
  using std::log;

  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::neg_binomial_2_log_rng(6, 2, rng));
  EXPECT_NO_THROW(stan::math::neg_binomial_2_log_rng(-0.5, 1, rng));
  EXPECT_NO_THROW(stan::math::neg_binomial_2_log_rng(log(1e8), 1, rng));

  EXPECT_THROW(stan::math::neg_binomial_2_log_rng(0, -2, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_rng(6, -2, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_rng(-6, -0.1, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_rng(
                   stan::math::positive_infinity(), 2, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_rng(
                   stan::math::positive_infinity(), 6, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_rng(
                   2, stan::math::positive_infinity(), rng),
               std::domain_error);

  std::string error_msg;

  error_msg
      = "neg_binomial_2_log_rng: Exponential "
        "of the log-location parameter divided by the precision "
        "parameter is inf";
  try {
    stan::math::neg_binomial_2_log_rng(log(1e300), 1e-300, rng);
    FAIL() << "neg_binomial_2_log_rng should have thrown" << std::endl;
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(error_msg) == std::string::npos)
      FAIL() << "Error message is different than expected" << std::endl
             << "EXPECTED: " << error_msg << std::endl
             << "FOUND: " << e.what() << std::endl;
    SUCCEED();
  }

  error_msg
      = "neg_binomial_2_log_rng: Random number that "
        "came from gamma distribution is";
  try {
    stan::math::neg_binomial_2_log_rng(log(1e10), 1e20, rng);
    FAIL() << "neg_binomial_2_log_rng should have thrown" << std::endl;
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(error_msg) == std::string::npos)
      FAIL() << "Error message is different than expected" << std::endl
             << "EXPECTED: " << error_msg << std::endl
             << "FOUND: " << e.what() << std::endl;
    SUCCEED();
  }
}

TEST(ProbDistributionsNegBinomial2Log, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::negative_binomial_distribution<> dist(1.1,
                                                     1.1 / (1.1 + exp(2.4)));
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
    int a = stan::math::neg_binomial_2_log_rng(2.4, 1.1, rng);
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

TEST(ProbDistributionsNegBinomial2Log, chiSquareGoodnessFitTest2) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::negative_binomial_distribution<> dist(0.6,
                                                     0.6 / (0.6 + exp(2.4)));
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
    int a = stan::math::neg_binomial_2_log_rng(2.4, 0.6, rng);
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

TEST(ProbDistributionsNegBinomial2Log, chiSquareGoodnessFitTest3) {
  boost::random::mt19937 rng;
  int N = 1000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::negative_binomial_distribution<> dist(121,
                                                     121 / (121 + exp(2.4)));
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
    int a = stan::math::neg_binomial_2_log_rng(2.4, 121, rng);
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

TEST(ProbNegBinomial2, log_matches_lpmf) {
  double y = 0.8;
  double mu = 1.1;
  double phi = 2.3;

  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_lpmf(y, mu, phi)),
                  (stan::math::neg_binomial_2_log(y, mu, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_lpmf<true>(y, mu, phi)),
                  (stan::math::neg_binomial_2_log<true>(y, mu, phi)));
  EXPECT_FLOAT_EQ((stan::math::neg_binomial_2_lpmf<false>(y, mu, phi)),
                  (stan::math::neg_binomial_2_log<false>(y, mu, phi)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_2_lpmf<true, double, double, double>(y, mu,
                                                                     phi)),
      (stan::math::neg_binomial_2_log<true, double, double, double>(y, mu,
                                                                    phi)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_2_lpmf<false, double, double, double>(y, mu,
                                                                      phi)),
      (stan::math::neg_binomial_2_log<false, double, double, double>(y, mu,
                                                                     phi)));
  EXPECT_FLOAT_EQ(
      (stan::math::neg_binomial_2_lpmf<double, double, double>(y, mu, phi)),
      (stan::math::neg_binomial_2_log<double, double, double>(y, mu, phi)));
}

TEST(ProbDistributionsNegBinomial2Log, neg_binomial_2_log_grid_test) {
  std::vector<double> mu_log_to_test
      = {-101, -27, -3, -1, -0.132, 0, 4, 10, 87};
  std::vector<double> phi_to_test = {2e-5, 0.36, 1, 10, 2.3e5, 1.8e10, 6e16};
  std::vector<int> n_to_test = {0, 1, 10, 39, 101, 3048, 150054};

  for (double mu_log : mu_log_to_test) {
    for (double phi : phi_to_test) {
      for (int n : n_to_test) {
        double val_log = stan::math::neg_binomial_2_log_lpmf(n, mu_log, phi);
        std::stringstream msg;
        double val_orig
            = stan::math::neg_binomial_2_lpmf(n, std::exp(mu_log), phi);
        msg << std::setprecision(22)
            << "neg_binomial_2_log_lpmf yields different result (" << val_log
            << ") than neg_binomial_2_lpmf (" << val_orig << ") for n = " << n
            << ", mu_log = " << mu_log << ", phi = " << phi << ".";
        stan::test::expect_near_rel(msg.str(), val_log, val_orig);
      }
    }
  }
}
