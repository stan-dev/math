#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, stdvec_corr_matrix_j) {
  size_t K = 4;
  size_t K_choose_2 = 6;
  Eigen::VectorXd x(K_choose_2);
  x << -1.0, 2.0, 0.0, 1.0, 3.0, -1.5;
  std::vector<Eigen::VectorXd> x_vec = {x, x, x};
  double lp = -12.9;
  std::vector<Eigen::MatrixXd> y_vec
      = stan::math::corr_matrix_constrain<true>(x_vec, K, lp);
  std::vector<Eigen::VectorXd> xrt_vec = stan::math::corr_matrix_free(y_vec);
  for (auto&& x_i : xrt_vec) {
    EXPECT_MATRIX_FLOAT_EQ(x, x_i);
  }
}

TEST(prob_transform, corr_matrix_j2x2) {
  // tests K=2 boundary case
  size_t K = 2;
  size_t K_choose_2 = 1;
  Eigen::VectorXd x(K_choose_2);
  x << -1.3;
  double lp = -12.9;
  Eigen::MatrixXd y = stan::math::corr_matrix_constrain(x, K, lp);
  Eigen::VectorXd xrt = stan::math::corr_matrix_free(y);
  EXPECT_MATRIX_FLOAT_EQ(x, xrt);
}

TEST(prob_transform, corr_matrix_j1x1) {
  // tests K=1, which has a different implementation
  size_t K = 1;
  size_t K_choose_2 = 0;
  Eigen::VectorXd x(K_choose_2);
  double lp = -12.9;
  Eigen::MatrixXd y = stan::math::corr_matrix_constrain(x, K, lp);
  Eigen::VectorXd xrt = stan::math::corr_matrix_free(y);
  EXPECT_MATRIX_FLOAT_EQ(x, xrt);
}

TEST(prob_transform, corr_matrix_constrain_exception) {
  unsigned int K = 4;
  unsigned int K_choose_2 = 6;
  Eigen::VectorXd x(K_choose_2 - 1);
  double lp = -12.9;

  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K), std::invalid_argument);
  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K, lp),
               std::invalid_argument);

  x.resize(K_choose_2 + 1);
  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K), std::invalid_argument);
  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K, lp),
               std::invalid_argument);
}

TEST(prob_transform, corr_matrix_rt) {
  unsigned int K = 4;
  unsigned int K_choose_2 = 6;
  Eigen::VectorXd x(K_choose_2);
  x << -1.0, 2.0, 0.0, 1.0, 3.0, -1.5;
  Eigen::MatrixXd y = stan::math::corr_matrix_constrain(x, K);
  Eigen::VectorXd xrt = stan::math::corr_matrix_free(y);
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}

TEST(prob_transform, corr_matrix_free_exception) {
  Eigen::MatrixXd y;

  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);
  y.resize(0, 10);
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);
  y.resize(10, 0);
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);
  y.resize(1, 2);
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);

  y.resize(2, 2);
  y << 0, 0, 0, 0;
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::domain_error);
}
