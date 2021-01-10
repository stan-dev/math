#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsCategoricalLogit, Categorical) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::log_softmax;
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << -1, 2, -10;
  Matrix<double, Dynamic, 1> theta_log_softmax = log_softmax(theta);

  EXPECT_FLOAT_EQ(theta_log_softmax[0],
                  stan::math::categorical_logit_log(1, theta));
  EXPECT_FLOAT_EQ(theta_log_softmax[1],
                  stan::math::categorical_logit_log(2, theta));
  EXPECT_FLOAT_EQ(theta_log_softmax[2],
                  stan::math::categorical_logit_log(3, theta));
}

TEST(ProbDistributionsCategoricalLogit, CategoricalVectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using stan::math::log_softmax;
  Matrix<double, Dynamic, 1> theta(3);
  theta << -1, 2, -10;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_log(ns, theta));

  Matrix<double, Dynamic, 1> theta_log_softmax = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(
      theta_log_softmax[0] + theta_log_softmax[1] + theta_log_softmax[0],
      stan::math::categorical_logit_log(ms, theta));

  std::vector<VectorXd> arr_theta(3);
  arr_theta[0] = log_softmax((VectorXd(3) << 1.5, -3, 8.2).finished());
  arr_theta[1] = log_softmax((VectorXd(3) << 5.1, 3.2, 5.2).finished());
  arr_theta[2] = log_softmax((VectorXd(3) << -3.1, -9.8, -2.1).finished());

  EXPECT_FLOAT_EQ(arr_theta[0][0] + arr_theta[1][1] + arr_theta[2][0],
                  stan::math::categorical_logit_log(ms, arr_theta));

  EXPECT_FLOAT_EQ(arr_theta[0][2] + arr_theta[1][2] + arr_theta[2][2],
                  stan::math::categorical_logit_log(3, arr_theta));
}

TEST(ProbDistributionsCategoricalLogit, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> theta(3, 1);
  theta << -1, 2, 10;
  EXPECT_FLOAT_EQ(0, stan::math::categorical_logit_log<true>(1, theta));
  EXPECT_FLOAT_EQ(0, stan::math::categorical_logit_log<true>(3, theta));
}

TEST(ProbDistributionsCategoricalLogit, error) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::categorical_logit_log;

  unsigned int n = 1;
  unsigned int N = 3;
  Matrix<double, Dynamic, 1> theta(N, 1);
  theta << 0.3, 0.5, 0.2;

  EXPECT_NO_THROW(categorical_logit_log(N, theta));
  EXPECT_NO_THROW(categorical_logit_log(n, theta));
  EXPECT_NO_THROW(categorical_logit_log(2, theta));
  EXPECT_THROW(categorical_logit_log(N + 1, theta), std::domain_error);
  EXPECT_THROW(categorical_logit_log(0, theta), std::domain_error);

  theta(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(categorical_logit_log(1, theta), std::domain_error);

  theta(1) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(categorical_logit_log(1, theta), std::domain_error);

  std::vector<int> ns(2);
  ns[0] = 1;
  ns[1] = 2;
  Eigen::VectorXd theta2(2);
  theta2 << 0.3, 0.5;
  EXPECT_NO_THROW(categorical_logit_log(ns, theta2));

  ns[0] = -1;
  EXPECT_THROW(categorical_logit_log(ns, theta2), std::domain_error);

  ns[0] = 1;
  ns[1] = 12;
  EXPECT_THROW(categorical_logit_log(ns, theta2), std::domain_error);
}

TEST(ProbDistributionsCategoricalLogit, error_vec) {
  Eigen::VectorXd beta1 = Eigen::VectorXd::Random(4);
  Eigen::VectorXd beta2 = Eigen::VectorXd::Random(4);
  std::vector<Eigen::VectorXd> betas = {beta1, beta2};

  double inf = std::numeric_limits<double>::infinity();

  // Check consistent sizes
  {
    std::vector<int> ns = {1};
    EXPECT_THROW(stan::math::categorical_logit_lpmf(ns, betas),
                 std::invalid_argument);
    ns = {1, 2, 3};
    EXPECT_THROW(stan::math::categorical_logit_lpmf(ns, betas),
                 std::invalid_argument);
    EXPECT_NO_THROW(stan::math::categorical_logit_lpmf(ns, beta1));
  }

  // Check bounded
  {
    std::vector<int> ns = {1, 10000};
    EXPECT_THROW(stan::math::categorical_logit_lpmf(ns, betas),
                 std::domain_error);
  }

  // Check finite
  {
    std::vector<Eigen::VectorXd> betas_inf = betas;
    betas_inf[1](0) = inf;
    EXPECT_THROW(stan::math::categorical_logit_lpmf(1, betas_inf),
                 std::domain_error);
  }
}
