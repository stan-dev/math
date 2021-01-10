#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>
#include <limits>

TEST(ProbDistributionsCategorical, Categorical) {
  using Eigen::Dynamic;
  using Eigen::VectorXd;
  VectorXd theta(3);
  theta << 0.3, 0.5, 0.2;
  EXPECT_FLOAT_EQ(-1.203973, stan::math::categorical_log(1, theta));
  EXPECT_FLOAT_EQ(-0.6931472, stan::math::categorical_log(2, theta));
}
TEST(ProbDistributionsCategorical, Propto) {
  using Eigen::Dynamic;
  using Eigen::VectorXd;
  VectorXd theta(3);
  theta << 0.3, 0.5, 0.2;
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_log<true>(1, theta));
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_log<true>(2, theta));
}

TEST(ProbDistributionsCategorical, VectorInt) {
  using Eigen::Dynamic;
  using Eigen::VectorXd;
  VectorXd theta(3);
  theta << 0.3, 0.5, 0.2;
  std::vector<int> xs0;
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_log(xs0, theta));

  std::vector<int> xs(3);
  xs[0] = 1;
  xs[1] = 3;
  xs[2] = 1;

  EXPECT_FLOAT_EQ(log(0.3) + log(0.2) + log(0.3),
                  stan::math::categorical_log(xs, theta));

  std::vector<VectorXd> arr_theta(3);
  arr_theta[0] = (VectorXd(3) << 0.1, 0.6, 0.3).finished();
  arr_theta[1] = (VectorXd(3) << 0.4, 0.2, 0.4).finished();
  arr_theta[2] = (VectorXd(3) << 0.5, 0.2, 0.3).finished();

  EXPECT_FLOAT_EQ(log(0.1) + log(0.4) + log(0.5),
                  stan::math::categorical_log(xs, arr_theta));

  EXPECT_FLOAT_EQ(log(0.6) + log(0.2) + log(0.2),
                  stan::math::categorical_log(2, arr_theta));
}

TEST(ProbDistributionsCategorical, error) {
  using Eigen::Dynamic;
  using Eigen::VectorXd;
  using stan::math::categorical_log;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  unsigned int n = 1;
  unsigned int N = 3;
  VectorXd theta(N);
  theta << 0.3, 0.5, 0.2;

  EXPECT_NO_THROW(categorical_log(N, theta));
  EXPECT_NO_THROW(categorical_log(n, theta));
  EXPECT_NO_THROW(categorical_log(2, theta));
  EXPECT_THROW(categorical_log(N + 1, theta), std::domain_error);
  EXPECT_THROW(categorical_log(0, theta), std::domain_error);

  theta(0) = nan;
  EXPECT_THROW(categorical_log(n, theta), std::domain_error);
  theta(0) = inf;
  EXPECT_THROW(categorical_log(n, theta), std::domain_error);
  theta(0) = -inf;
  EXPECT_THROW(categorical_log(n, theta), std::domain_error);
  theta(0) = -1;
  theta(1) = 1;
  theta(2) = 0;
  EXPECT_THROW(categorical_log(n, theta), std::domain_error);

  std::vector<int> ns(3);
  ns[0] = 3;
  ns[1] = 2;
  ns[2] = 3;
  EXPECT_THROW(categorical_log(ns, theta), std::domain_error);

  theta << 0.3, 0.5, 0.2;
  EXPECT_NO_THROW(categorical_log(ns, theta));
  ns[1] = -1;
  EXPECT_THROW(categorical_log(ns, theta), std::domain_error);

  ns[1] = 1;
  ns[2] = 12;
  EXPECT_THROW(categorical_log(ns, theta), std::domain_error);
}

TEST(ProbDistributionsCategorical, error_check) {
  using Eigen::Dynamic;
  using Eigen::VectorXd;
  using stan::math::categorical_log;
  boost::random::mt19937 rng;

  VectorXd theta(3);
  theta << 0.15, 0.45, 0.50;

  EXPECT_THROW(stan::math::categorical_rng(theta, rng), std::domain_error);
}

TEST(ProbDistributionsCategorical, chiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::VectorXd;
  using stan::math::categorical_log;
  boost::random::mt19937 rng;

  int N = 10000;
  VectorXd theta(3);
  theta << 0.15, 0.45, 0.40;
  int K = theta.rows();
  boost::math::chi_squared mydist(K - 1);

  VectorXd loc(theta.rows());
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

TEST(ProbDistributionsCategorical, categorical_vecRNG) {
  using stan::math::categorical_rng;
  using stan::math::softmax;
  boost::random::mt19937 rng;

  Eigen::VectorXd vec1(4);
  vec1 << softmax(Eigen::VectorXd::Random(4));
  Eigen::VectorXd vec2(4);
  vec2 << softmax(Eigen::VectorXd::Random(4));
  Eigen::VectorXd vec3(4);
  vec3 << softmax(Eigen::VectorXd::Random(4));

  std::vector<Eigen::VectorXd> vecs{vec1, vec2, vec3};

  std::vector<int> rng_stvec = categorical_rng(vecs, rng);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(rng_stvec[i], categorical_rng(vecs[i], rng));
  }
}

TEST(ProbDistributionsCategorical, categorical_vecRNG_throw) {
  using stan::math::categorical_rng;
  using stan::math::softmax;
  boost::random::mt19937 rng;

  Eigen::VectorXd vec1(4);
  vec1 << softmax(Eigen::VectorXd::Random(4));
  Eigen::VectorXd vec2 = Eigen::VectorXd::Zero(4);
  vec2(0) = 2.0;
  Eigen::VectorXd vec3(4);
  vec3 << softmax(Eigen::VectorXd::Random(4));

  std::vector<Eigen::VectorXd> vecs{vec1, vec2, vec3};

  EXPECT_THROW(categorical_rng(vecs, rng), std::domain_error);
}

TEST(ProbDistributionsCategorical, error_vec) {
  Eigen::VectorXd beta1 = stan::math::softmax(Eigen::VectorXd::Random(4));
  Eigen::VectorXd beta2 = stan::math::softmax(Eigen::VectorXd::Random(4));
  std::vector<Eigen::VectorXd> betas = {beta1, beta2};

  double inf = std::numeric_limits<double>::infinity();

  // Check consistent sizes
  {
    std::vector<int> ns = {1};
    EXPECT_THROW(stan::math::categorical_lpmf(ns, betas),
                 std::invalid_argument);
    ns = {1, 2, 3};
    EXPECT_THROW(stan::math::categorical_lpmf(ns, betas),
                 std::invalid_argument);
    EXPECT_NO_THROW(stan::math::categorical_lpmf(ns, beta1));
  }

  // Check bounded
  {
    std::vector<int> ns = {1, 10000};
    EXPECT_THROW(stan::math::categorical_lpmf(ns, betas), std::domain_error);
  }

  // Check simplex
  {
    std::vector<Eigen::VectorXd> betas_not_simplex = betas;
    betas_not_simplex[1](0) = 2.0;
    betas_not_simplex[1](1) = 0.0;
    betas_not_simplex[1](2) = 0.0;
    EXPECT_THROW(stan::math::categorical_lpmf(1, betas_not_simplex),
                 std::domain_error);
  }
}
