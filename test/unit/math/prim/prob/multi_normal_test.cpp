#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsMultiNormal, NotVectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;

  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  EXPECT_FLOAT_EQ(-11.73908, stan::math::multi_normal_lpdf(y, mu, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(mu, Sigma, rng));
}

TEST(ProbDistributionsMultiNormal, Vectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  vector<Matrix<double, Dynamic, 1> > vec_y(2);
  vector<Matrix<double, 1, Dynamic> > vec_y_t(2);
  Matrix<double, Dynamic, 1> y(3);
  Matrix<double, 1, Dynamic> y_t(3);
  y << 2.0, -2.0, 11.0;
  vec_y[0] = y;
  vec_y_t[0] = y;
  y << 4.0, -2.0, 1.0;
  vec_y[1] = y;
  vec_y_t[1] = y;
  y_t = y;

  vector<Matrix<double, Dynamic, 1> > vec_mu(2);
  vector<Matrix<double, 1, Dynamic> > vec_mu_t(2);
  Matrix<double, Dynamic, 1> mu(3);
  Matrix<double, 1, Dynamic> mu_t(3);
  mu << 1.0, -1.0, 3.0;
  vec_mu[0] = mu;
  vec_mu_t[0] = mu;
  mu << 2.0, -1.0, 4.0;
  vec_mu[1] = mu;
  vec_mu_t[1] = mu;
  mu_t = mu;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 10.0, -3.0, 0.0, -3.0, 5.0, 0.0, 0.0, 0.0, 5.0;

  // y and mu vectorized
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_lpdf(vec_y, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_lpdf(vec_y_t, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_lpdf(vec_y, vec_mu_t, Sigma));
  EXPECT_FLOAT_EQ(-11.928077 - 6.5378327,
                  stan::math::multi_normal_lpdf(vec_y_t, vec_mu_t, Sigma));

  // y vectorized
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_lpdf(vec_y, mu, Sigma));
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_lpdf(vec_y_t, mu, Sigma));
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_lpdf(vec_y, mu_t, Sigma));
  EXPECT_FLOAT_EQ(-10.44027 - 6.537833,
                  stan::math::multi_normal_lpdf(vec_y_t, mu_t, Sigma));

  // mu vectorized
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_lpdf(y, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_lpdf(y_t, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_lpdf(y, vec_mu_t, Sigma));
  EXPECT_FLOAT_EQ(-6.26954 - 6.537833,
                  stan::math::multi_normal_lpdf(y_t, vec_mu_t, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(vec_mu, Sigma, rng));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(vec_mu_t, Sigma, rng));
}
TEST(ProbDistributionsMultiNormal, Sigma) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(2, 1);
  y << 2.0, -2.0;
  Matrix<double, Dynamic, 1> mu(2, 1);
  mu << 1.0, -1.0;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 9.0, -3.0, -3.0, 4.0;
  EXPECT_NO_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(mu, Sigma, rng));

  // non-symmetric
  Sigma(0, 1) = -2.5;
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
}
TEST(ProbDistributionsMultiNormal, Mu) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  EXPECT_NO_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(mu, Sigma, rng));

  mu(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
  mu(0) = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
  mu(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
}
TEST(ProbDistributionsMultiNormal, MultiNormalOneRow) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, 1, Dynamic> y(1, 3);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  EXPECT_FLOAT_EQ(-11.73908, stan::math::multi_normal_lpdf(y, mu, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(mu, Sigma, rng));
}

TEST(ProbDistributionsMultiNormal, SigmaMultiRow) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, 1, Dynamic> y(1, 2);
  y << 2.0, -2.0;
  Matrix<double, Dynamic, 1> mu(2, 1);
  mu << 1.0, -1.0;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 9.0, -3.0, -3.0, 4.0;
  EXPECT_NO_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(mu, Sigma, rng));

  // wrong dimensions
  Matrix<double, Dynamic, 1> z(3, 1);
  z << 2.0, -2.0, 1.0;
  EXPECT_THROW(stan::math::multi_normal_lpdf(z, mu, Sigma),
               std::invalid_argument);

  // non-symmetric
  Sigma(0, 1) = -2.5;
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
}
TEST(ProbDistributionsMultiNormal, MuMultiRow) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, 1, Dynamic> y(1, 3);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  EXPECT_NO_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_rng(mu, Sigma, rng));

  mu(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
  mu(0) = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
  mu(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma), std::domain_error);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng), std::domain_error);
}
TEST(ProbDistributionsMultiNormal, SizeMismatch) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, 1, Dynamic> y(1, 3);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(2, 1);
  mu << 1.0, -1.0;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0;
  EXPECT_THROW(stan::math::multi_normal_lpdf(y, mu, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_rng(mu, Sigma, rng),
               std::invalid_argument);
}

TEST(ProbDistributionsMultiNormal, marginalOneChiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, Dynamic> sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 3.0;
  std::vector<Matrix<double, Dynamic, 1> > mu(3);
  mu[0].resize(3);
  mu[1].resize(3);
  mu[2].resize(3);
  mu[0] << 2.0, -2.0, 11.0;
  mu[1] << 7.0, -3.0, 5.0;
  mu[2] << 5.0, -6.0, 1.0;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(2.0, 3.0);
  boost::math::chi_squared mydist(K - 1);

  double loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = quantile(dist, i * std::pow(K, -1.0));

  int count = 0;
  int bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N / K;
  }
  Eigen::VectorXd a(mu[0].rows());
  while (count < N) {
    a = stan::math::multi_normal_rng(mu, sigma, rng)[0];
    int i = 0;
    while (i < K - 1 && a(0) > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsMultiNormal, marginalTwoChiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, Dynamic> sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 3.0;
  std::vector<Matrix<double, 1, Dynamic> > mu(3);
  mu[0].resize(3);
  mu[1].resize(3);
  mu[2].resize(3);
  mu[0] << 1.0, 5.0, -1.0;
  mu[1] << 2.0, -2.0, 11.0;
  mu[2] << 7.0, 0.0, 3.0;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(-2.0, 2.0);
  boost::math::chi_squared mydist(K - 1);

  double loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = quantile(dist, i * std::pow(K, -1.0));

  int count = 0;
  int bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N / K;
  }
  Eigen::VectorXd a(mu[0].rows());
  while (count < N) {
    a = stan::math::multi_normal_rng(mu, sigma, rng)[1];
    int i = 0;
    while (i < K - 1 && a(1) > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsMultiNormal, marginalThreeChiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, Dynamic> sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 16.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 2.0, -2.0, 11.0;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(11.0, 4.0);
  boost::math::chi_squared mydist(K - 1);

  double loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = quantile(dist, i * std::pow(K, -1.0));

  int count = 0;
  int bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N / K;
  }
  Eigen::VectorXd a(mu.rows());
  while (count < N) {
    a = stan::math::multi_normal_rng(mu, sigma, rng);
    int i = 0;
    while (i < K - 1 && a(2) > loc[i])
      ++i;
    ++bin[i];
    count++;
  }

  double chi = 0;
  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(multiNormalRng, nonPosDefErrorTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_normal_rng;
  using std::vector;
  Eigen::MatrixXd S(2, 2);
  S << 0, 1, 1, 0;  // not pos definite
  Eigen::VectorXd mu(2);
  mu << 1, 2;
  boost::random::mt19937 rng;
  EXPECT_THROW(multi_normal_rng(mu, S, rng), std::domain_error);
}

TEST(ProbDistributionsMultiNormal, WrongSize) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  vector<Matrix<double, Dynamic, 1> > y_3_3(3);
  vector<Matrix<double, Dynamic, 1> > y_3_1(3);
  vector<Matrix<double, Dynamic, 1> > y_3_2(3);
  vector<Matrix<double, Dynamic, 1> > y_1_3(1);
  vector<Matrix<double, Dynamic, 1> > y_2_3(2);
  Matrix<double, Dynamic, 1> y_3(3);
  Matrix<double, Dynamic, 1> y_2(2);
  Matrix<double, Dynamic, 1> y_1(1);
  y_3 << 2.0, -2.0, 11.0;
  y_2 << 2.0, -2.0;
  y_1 << 2.0;
  y_3_3[0] = y_3;
  y_3_3[1] = y_3;
  y_3_3[2] = y_3;
  y_3_1[0] = y_1;
  y_3_1[1] = y_1;
  y_3_1[2] = y_1;
  y_3_2[0] = y_2;
  y_3_2[1] = y_2;
  y_3_2[2] = y_2;
  y_1_3[0] = y_3;
  y_2_3[0] = y_3;
  y_2_3[1] = y_3;

  vector<Matrix<double, Dynamic, 1> > mu_3_3(3);
  Matrix<double, Dynamic, 1> mu_3(3);
  mu_3 << 2.0, -2.0, 11.0;
  mu_3_3[0] = mu_3;
  mu_3_3[1] = mu_3;
  mu_3_3[2] = mu_3;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 10.0, -3.0, 0.0, -3.0, 5.0, 0.0, 0.0, 0.0, 5.0;

  EXPECT_NO_THROW(stan::math::multi_normal_lpdf(y_3_3, mu_3_3, Sigma));
  EXPECT_NO_THROW(stan::math::multi_normal_lpdf(y_3, mu_3_3, Sigma));

  EXPECT_THROW(stan::math::multi_normal_lpdf(y_1_3, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_lpdf(y_2_3, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_lpdf(y_3_1, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_normal_lpdf(y_3_2, mu_3_3, Sigma),
               std::invalid_argument);
}
