#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>

stan::math::vector_d get_simplex_Phi(double lambda,
                                     const stan::math::vector_d& c) {
  using stan::math::Phi;
  int K = c.size() + 1;
  stan::math::vector_d theta(K);
  theta(0) = 1.0 - Phi(lambda - c(0));
  for (int k = 1; k < (K - 1); ++k)
    theta(k) = Phi(lambda - c(k - 1)) - Phi(lambda - c(k));
  // - 0.0
  theta(K - 1) = Phi(lambda - c(K - 2));
  return theta;
}

TEST(ProbDistributions, ordered_probit_stability) {
  using stan::math::is_inf;
  using stan::math::ordered_probit_log;

  Eigen::VectorXd c(3);
  c << -0.3, 0.1, 1.2;

  EXPECT_FALSE(is_inf(ordered_probit_log(1, 10, c)));
  EXPECT_FALSE(is_inf(ordered_probit_log(2, 10, c)));
  EXPECT_NE(ordered_probit_log(4, 10, c), 0);

  EXPECT_NE(ordered_probit_log(1, -38, c), 0);
  EXPECT_FALSE(is_inf(ordered_probit_log(2, -38, c)));
  EXPECT_FALSE(is_inf(ordered_probit_log(4, -38, c)));
}

TEST(ProbDistributions, ordered_probit_vals) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using stan::math::ordered_probit_log;
  using stan::math::Phi;

  int K = 5;
  Matrix<double, Dynamic, 1> c(K - 1);
  c << -1.7, -0.3, 1.2, 2.6;
  double lambda = 1.1;

  stan::math::vector_d theta = get_simplex_Phi(lambda, c);

  double sum = 0.0;
  for (int k = 0; k < theta.size(); ++k)
    sum += theta(k);
  EXPECT_FLOAT_EQ(1.0, sum);

  for (int k = 0; k < K; ++k)
    EXPECT_FLOAT_EQ(log(theta(k)), ordered_probit_log(k + 1, lambda, c));

  EXPECT_THROW(ordered_probit_log(0, lambda, c), std::domain_error);
  EXPECT_THROW(ordered_probit_log(6, lambda, c), std::domain_error);
}

TEST(ProbDistributions, ordered_probit_vals_2) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  using stan::math::ordered_probit_log;
  using stan::math::Phi;

  int K = 3;
  Matrix<double, Dynamic, 1> c(K - 1);
  c << -0.2, 4;
  double lambda = -0.9;

  stan::math::vector_d theta = get_simplex_Phi(lambda, c);

  double sum = 0.0;
  for (int k = 0; k < theta.size(); ++k)
    sum += theta(k);
  EXPECT_FLOAT_EQ(1.0, sum);

  for (int k = 0; k < K; ++k)
    EXPECT_FLOAT_EQ(log(theta(k)), ordered_probit_log(k + 1, lambda, c));

  EXPECT_THROW(ordered_probit_log(0, lambda, c), std::domain_error);
  EXPECT_THROW(ordered_probit_log(4, lambda, c), std::domain_error);
}

TEST(ProbDistributions, ordered_probit) {
  using stan::math::ordered_probit_log;
  int K = 4;
  Eigen::Matrix<double, Eigen::Dynamic, 1> c(K - 1);
  c << -0.3, 0.1, 1.2;
  double lambda = 0.5;
  EXPECT_THROW(ordered_probit_log(-1, lambda, c), std::domain_error);
  EXPECT_THROW(ordered_probit_log(0, lambda, c), std::domain_error);
  EXPECT_THROW(ordered_probit_log(5, lambda, c), std::domain_error);
  for (int k = 1; k <= K; ++k)
    EXPECT_NO_THROW(ordered_probit_log(k, lambda, c));

  // init size zero
  Eigen::Matrix<double, Eigen::Dynamic, 1> c_zero;
  EXPECT_EQ(0, c_zero.size());
  EXPECT_THROW(ordered_probit_log(1, lambda, c_zero), std::invalid_argument);

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_neg(1);
  c_neg << -13.7;
  EXPECT_NO_THROW(ordered_probit_log(1, lambda, c_neg));

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_unord(3);
  c_unord << 1.0, 0.4, 2.0;
  EXPECT_THROW(ordered_probit_log(1, lambda, c_unord), std::domain_error);

  Eigen::Matrix<double, Eigen::Dynamic, 1> c_unord_2(3);
  c_unord_2 << 1.0, 2.0, 0.4;
  EXPECT_THROW(ordered_probit_log(1, lambda, c_unord_2), std::domain_error);

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  EXPECT_THROW(ordered_probit_log(1, nan, c), std::domain_error);
  EXPECT_THROW(ordered_probit_log(1, inf, c), std::domain_error);

  Eigen::Matrix<double, Eigen::Dynamic, 1> cbad(2);
  cbad << 0.2, inf;
  EXPECT_THROW(ordered_probit_log(1, 1.0, cbad), std::domain_error);
  cbad[1] = nan;
  EXPECT_THROW(ordered_probit_log(1, 1.0, cbad), std::domain_error);

  Eigen::Matrix<double, Eigen::Dynamic, 1> cbad1(1);
  cbad1 << inf;
  EXPECT_THROW(ordered_probit_log(1, 1.0, cbad1), std::domain_error);
  cbad1[0] = nan;
  EXPECT_THROW(ordered_probit_log(1, 1.0, cbad1), std::domain_error);

  Eigen::Matrix<double, Eigen::Dynamic, 1> cbad3(3);
  cbad3 << 0.5, inf, 1.0;
  EXPECT_THROW(ordered_probit_log(1, 1.0, cbad3), std::domain_error);
  cbad3[1] = nan;
  EXPECT_THROW(ordered_probit_log(1, 1.0, cbad3), std::domain_error);
}

TEST(ProbDistributionOrderedProbit, error_check) {
  boost::random::mt19937 rng;
  double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd c(4);
  c << -2, 2.0, 5, 10;
  EXPECT_NO_THROW(stan::math::ordered_probit_rng(4.0, c, rng));

  EXPECT_THROW(
      stan::math::ordered_probit_rng(stan::math::positive_infinity(), c, rng),
      std::domain_error);
  c << -inf, 2.0, -5, inf;
  EXPECT_THROW(stan::math::ordered_probit_rng(4.0, c, rng), std::domain_error);
}

TEST(ProbDistributionOrderedProbit, chiSquareGoodnessFitTest) {
  using stan::math::Phi;
  boost::random::mt19937 rng;
  int N = 10000;
  double eta = 1.0;
  Eigen::VectorXd theta(3);
  theta << -0.4, 4.0, 6.2;
  Eigen::VectorXd prob(4);
  prob(0) = 1 - Phi(eta - theta(0));
  prob(1) = Phi(eta - theta(0)) - Phi(eta - theta(1));
  prob(2) = Phi(eta - theta(1)) - Phi(eta - theta(2));
  prob(3) = Phi(eta - theta(2));
  int K = prob.rows();
  boost::math::chi_squared mydist(K - 1);

  Eigen::VectorXd loc(prob.rows());
  for (int i = 0; i < prob.rows(); i++)
    loc(i) = 0;

  for (int i = 0; i < prob.rows(); i++) {
    for (int j = i; j < prob.rows(); j++)
      loc(j) += prob(i);
  }

  int count = 0;
  int bin[K];
  double expect[K];
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = N * prob(i);
  }

  while (count < N) {
    int a = stan::math::ordered_probit_rng(eta, theta, rng);
    bin[a - 1]++;
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++)
    chi += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);
  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
