#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>
#include <limits>

TEST(ProbDistributionsCategorical, Categorical) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 0.3, 0.5, 0.2;
  EXPECT_FLOAT_EQ(-1.203973, stan::math::categorical_lpmf(1, theta));
  EXPECT_FLOAT_EQ(-0.6931472, stan::math::categorical_lpmf(2, theta));
}
TEST(ProbDistributionsCategorical, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 0.3, 0.5, 0.2;
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_lpmf<true>(1, theta));
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_lpmf<true>(2, theta));
}

TEST(ProbDistributionsCategorical, VectorInt) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << 0.3, 0.5, 0.2;
  std::vector<int> xs0;
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_lpmf(xs0, theta));

  std::vector<int> xs(3);
  xs[0] = 1;
  xs[1] = 3;
  xs[2] = 1;

  EXPECT_FLOAT_EQ(log(0.3) + log(0.2) + log(0.3),
                  stan::math::categorical_lpmf(xs, theta));
}

TEST(ProbDistributionsCategorical, error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::categorical_lpmf;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  unsigned int n = 1;
  unsigned int N = 3;
  Matrix<double, Dynamic, 1> theta(N, 1);
  theta << 0.3, 0.5, 0.2;

  EXPECT_NO_THROW(categorical_lpmf(N, theta));
  EXPECT_NO_THROW(categorical_lpmf(n, theta));
  EXPECT_NO_THROW(categorical_lpmf(2, theta));
  EXPECT_THROW(categorical_lpmf(N + 1, theta), std::domain_error);
  EXPECT_THROW(categorical_lpmf(0, theta), std::domain_error);

  theta(0) = nan;
  EXPECT_THROW(categorical_lpmf(n, theta), std::domain_error);
  theta(0) = inf;
  EXPECT_THROW(categorical_lpmf(n, theta), std::domain_error);
  theta(0) = -inf;
  EXPECT_THROW(categorical_lpmf(n, theta), std::domain_error);
  theta(0) = -1;
  theta(1) = 1;
  theta(2) = 0;
  EXPECT_THROW(categorical_lpmf(n, theta), std::domain_error);

  std::vector<int> ns(3);
  ns[0] = 3;
  ns[1] = 2;
  ns[2] = 3;
  EXPECT_THROW(categorical_lpmf(ns, theta), std::domain_error);

  theta << 0.3, 0.5, 0.2;
  EXPECT_NO_THROW(categorical_lpmf(ns, theta));
  ns[1] = -1;
  EXPECT_THROW(categorical_lpmf(ns, theta), std::domain_error);

  ns[1] = 1;
  ns[2] = 12;
  EXPECT_THROW(categorical_lpmf(ns, theta), std::domain_error);
}

TEST(ProbDistributionsCategorical, error_check) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::categorical_lpmf;
  boost::random::mt19937 rng;

  Matrix<double, Dynamic, Dynamic> theta(3, 1);
  theta << 0.15, 0.45, 0.50;

  EXPECT_THROW(stan::math::categorical_rng(theta, rng), std::domain_error);
}

TEST(ProbDistributionsCategorical, chiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::categorical_lpmf;
  boost::random::mt19937 rng;

  int N = 10000;
  Matrix<double, Dynamic, Dynamic> theta(3, 1);
  theta << 0.15, 0.45, 0.40;
  int K = theta.rows();
  boost::math::chi_squared mydist(K - 1);

  Eigen::Matrix<double, Eigen::Dynamic, 1> loc(theta.rows(), 1);
  for (int i = 0; i < theta.rows(); i++)
    loc(i) = 0;

  for (int i = 0; i < theta.rows(); i++) {
    for (int j = i; j < theta.rows(); j++)
      loc(j) += theta(i);
  }

  int count = 0;
  int bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N * theta(i);
  }

  while (count < N) {
    int a = stan::math::categorical_rng(theta, rng);
    bin[a - 1]++;
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
