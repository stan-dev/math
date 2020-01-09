#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

class BetaProportionTestRig : public VectorRealRNGTestRig {
 public:
  BetaProportionTestRig()
      : VectorRealRNGTestRig(10000, 10, {0.3, 0.4, 0.5, 0.6, 0.7}, {1, 2, 3},
                             {-2.5, -1.7, -0.1, 0.0}, {-3, -2, -1, 0},
                             {0.35, 0.5, 0.9, 1.7, 2.1, 4.1}, {1, 2, 3, 4},
                             {-2.7, -1.5, -0.5, 0.0}, {-3, -2, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& kappa, const T3&,
                        T_rng& rng) const {
    return stan::math::beta_proportion_rng(mu, kappa, rng);
  }

  std::vector<double> generate_quantiles(double mu, double kappa,
                                         double) const {
    // transform from location and precision parameterization
    // into shape1 (alpha) and shape2 (beta) parameterization
    double alpha = mu * kappa;
    double beta = kappa - alpha;
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::beta_distribution<> dist(alpha, beta);

    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());

    return quantiles;
  }
};

TEST(ProbDistributionsBetaProportion, errorCheck) {
  check_dist_throws_real_first_argument(BetaProportionTestRig());
}

TEST(ProbDistributionsBetaProportion, distributionTest) {
  check_quantiles_real_first_argument(BetaProportionTestRig());
}

TEST(ProbDistributionsBetaProportion, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::beta_proportion_rng(0.5, 3.0, rng));
  EXPECT_NO_THROW(stan::math::beta_proportion_rng(1e-10, 1e-10, rng));
  EXPECT_THROW(stan::math::beta_proportion_rng(2.0, -1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_rng(-2.0, 1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_rng(-2.0, -1.0, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_rng(stan::math::positive_infinity(),
                                               1.0, rng),
               std::domain_error);
  EXPECT_THROW(
      stan::math::beta_proportion_rng(2, stan::math::positive_infinity(), rng),
      std::domain_error);
}

TEST(ProbDistributionsBetaProportion, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  double mu = 0.5;     // location
  double kappa = 3.0;  // precision

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::beta_proportion_rng(mu, kappa, rng));
  }

  // transform from location and precision parameterization
  // into shape1 (alpha) and shape2 (beta) parameterization
  double alpha = mu * kappa;
  double beta = kappa - alpha;

  // Generate quantiles from boost's beta distribution
  boost::math::beta_distribution<> dist(alpha, beta);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}

TEST(ProbDistributionsBetaProportion, chiSquareGoodnessFitTest2) {
  boost::random::mt19937 rng;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));

  double mu = 0.3;     // location
  double kappa = 0.5;  // precision

  std::vector<double> samples;
  for (int i = 0; i < N; ++i) {
    samples.push_back(stan::math::beta_proportion_rng(mu, kappa, rng));
  }

  // transform from location and precision parameterization
  // into shape1 (alpha) and shape2 (beta) parameterization
  double alpha = mu * kappa;
  double beta = kappa - alpha;

  // Generate quantiles from boost's beta distribution
  boost::math::beta_distribution<> dist(alpha, beta);
  std::vector<double> quantiles;
  for (int i = 1; i < K; ++i) {
    double frac = static_cast<double>(i) / K;
    quantiles.push_back(quantile(dist, frac));
  }
  quantiles.push_back(std::numeric_limits<double>::max());

  // Assert that they match
  assert_matches_quantiles(samples, quantiles, 1e-6);
}
