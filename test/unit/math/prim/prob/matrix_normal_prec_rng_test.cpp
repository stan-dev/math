#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <exception>
#include <limits>
#include <vector>

TEST(ProbDistributionsMatrixNormalPrecRng, ErrorSigma) {
  using Eigen::MatrixXd;
  using stan::math::matrix_normal_prec_rng;
  boost::random::mt19937 rng;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -std::numeric_limits<double>::infinity();

  MatrixXd Mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  MatrixXd D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  EXPECT_NO_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng));

  // non-symmetric
  Sigma(0, 1) = -2.5;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  Sigma(0, 1) = Sigma(1, 0);

  // non-spd
  Sigma(0, 0) = -3.0;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  Sigma(0, 0) = 9.0;

  // NaN
  Sigma(0, 0) = nan;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  Sigma(0, 0) = 9.0;

  // inf
  Sigma(0, 0) = inf;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  Sigma(0, 0) = ninf;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  Sigma(0, 0) = 9.0;

  // non square
  MatrixXd rect_Sigma(2, 3);
  rect_Sigma << 1, 0, 0, 0, 1, 0;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, rect_Sigma, D, rng),
               std::invalid_argument);
}

TEST(ProbDistributionsMatrixNormalPrecRng, ErrorD) {
  using Eigen::MatrixXd;
  using stan::math::matrix_normal_prec_rng;
  boost::random::mt19937 rng;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  double ninf = -std::numeric_limits<double>::infinity();

  MatrixXd Mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  MatrixXd D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  EXPECT_NO_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng));

  // non-symmetric
  D(0, 1) = -2.5;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  D(0, 1) = Sigma(1, 0);

  // non-spd
  D(0, 0) = -3.0;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  D(0, 0) = 1.0;

  // NaN
  D(0, 0) = nan;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  D(0, 0) = 1.0;

  // inf
  D(0, 0) = inf;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  D(0, 0) = ninf;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, D, rng), std::domain_error);
  D(0, 0) = 1.0;

  // non square
  MatrixXd rect_D(2, 3);
  rect_D << 1, 0, 0, 0, 1, 0;
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, rect_D, rng),
               std::invalid_argument);
}

TEST(ProbDistributionsMatrixNormalPrecRng, ErrorSize) {
  using Eigen::MatrixXd;
  using stan::math::matrix_normal_prec_rng;
  boost::random::mt19937 rng;

  MatrixXd Mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  MatrixXd SigmaWrong(4, 4);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;
  SigmaWrong << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  MatrixXd D(5, 5);
  MatrixXd DWrong(4, 4);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;
  DWrong << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

  EXPECT_THROW(matrix_normal_prec_rng(Mu, SigmaWrong, D, rng),
               std::invalid_argument);
  EXPECT_THROW(matrix_normal_prec_rng(Mu, Sigma, DWrong, rng),
               std::invalid_argument);
}

/**
 * Assert that the samples come from the normal distribution with this
 * mean and variance.
 */
void assert_matches_normal_distribution(const double mean,
                                        const double variance,
                                        const std::vector<double> &samples) {
  using Eigen::MatrixXd;
  using stan::math::matrix_normal_prec_rng;
  int N = samples.size();
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(mean, sqrt(variance));
  std::vector<double> quantiles;
  for (int i = 1; i < K; i++)
    quantiles.push_back(quantile(dist, i * std::pow(K, -1.0)));
  quantiles.push_back(std::numeric_limits<double>::max());
  assert_matches_quantiles(samples, quantiles, 1e-6);
}

std::vector<double> extract_entry(const unsigned int r, const unsigned int c,
                                  const std::vector<Eigen::MatrixXd> &samples) {
  using Eigen::MatrixXd;
  using stan::math::matrix_normal_prec_rng;
  std::vector<double> univariate_samples;
  for (auto sample : samples)
    univariate_samples.push_back(sample(r, c));
  return univariate_samples;
}

std::vector<double> extract_sum_of_entries(
    const unsigned int r1, const unsigned int c1, const unsigned int r2,
    const unsigned int c2, const std::vector<Eigen::MatrixXd> &samples) {
  using Eigen::MatrixXd;
  using stan::math::matrix_normal_prec_rng;
  std::vector<double> univariate_samples;
  for (auto sample : samples)
    univariate_samples.push_back(sample(r1, c1) + sample(r2, c2));
  return univariate_samples;
}

TEST(ProbDistributionsMatrixNormalPrecRng, marginalChiSquareGoodnessFitTest) {
  using Eigen::MatrixXd;
  using stan::math::matrix_normal_prec_rng;
  boost::random::mt19937 rng;

  MatrixXd Mu(2, 3);
  Mu << 1, 2, 3, 4, 5, 6;

  MatrixXd Sigma_cov(2, 2);
  Sigma_cov << 2, 1, 1, 3;
  MatrixXd Sigma = Sigma_cov.inverse();

  MatrixXd D_cov(3, 3);
  D_cov << 4, 1, 2, 1, 5, 3, 2, 3, 6;
  MatrixXd D = D_cov.inverse();

  // Kronecker product D_cov times Sigma_cov.
  // sage: Sigma_cov = Matrix([[2, 1], [1, 3]])
  // sage: D_cov = Matrix([[4, 1, 2], [1, 5, 3], [2, 3, 6]])
  // sage: D_Sigma_cov = D_cov.tensor_product(Sigma_cov)
  // sage: D_Sigma_cov
  //
  // [ 8  4| 2  1| 4  2]
  // [ 4 12| 1  3| 2  6]
  // [-----+-----+-----]
  // [ 2  1|10  5| 6  3]
  // [ 1  3| 5 15| 3  9]
  // [-----+-----+-----]
  // [ 4  2| 6  3|12  6]
  // [ 2  6| 3  9| 6 18]
  //
  // flatten(Mu) = [1, 4, 2, 5, 3, 6] ("vectorized".)
  // flatten(Y) ~ normal(flatten(Mu), D_Sigma_cov)
  // Where normal is the ordinary multivariate normal.

  int N = 10000;

  std::vector<Eigen::MatrixXd> samples;
  for (int count = 0; count < N; ++count) {
    MatrixXd Y = matrix_normal_prec_rng(Mu, Sigma, D, rng);
    samples.push_back(Y);
  }

  assert_matches_normal_distribution(1, 8, extract_entry(0, 0, samples));
  assert_matches_normal_distribution(2, 10, extract_entry(0, 1, samples));
  assert_matches_normal_distribution(3, 12, extract_entry(0, 2, samples));
  assert_matches_normal_distribution(4, 12, extract_entry(1, 0, samples));
  assert_matches_normal_distribution(5, 15, extract_entry(1, 1, samples));
  assert_matches_normal_distribution(6, 18, extract_entry(1, 2, samples));

  assert_matches_normal_distribution(
      // E[X + X2] = E[X] + E[X2].
      1 + 6,
      // var[X + X2] = var[X] + var[X2] + 2cov[X, X2].
      8 + 18 + 2 * 2, extract_sum_of_entries(0, 0, 1, 2, samples));
  assert_matches_normal_distribution(
      2 + 4, 10 + 12 + 2 * 1, extract_sum_of_entries(0, 1, 1, 0, samples));
}
