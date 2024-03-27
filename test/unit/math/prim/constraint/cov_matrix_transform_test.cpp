#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, cov_matrix_rt) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  unsigned int K = 4;
  double lp = 0;
  unsigned int K_choose_2 = 6;
  Matrix<double, Dynamic, 1> x(K_choose_2 + K);
  x << -1.0, 2.0, 0.0, 1.0, 3.0, -1.5, 1.0, 2.0, -1.5, 2.5;
  std::vector<Matrix<double, Dynamic, 1>> x_vec = {x, x, x};
  std::vector<Matrix<double, Dynamic, Dynamic>> y_vec
      = stan::math::cov_matrix_constrain<false>(x_vec, K, lp);
  std::vector<Eigen::VectorXd> xrt = stan::math::cov_matrix_free(y_vec);
  for (auto&& x_i : xrt) {
    EXPECT_EQ(x.size(), x_i.size());
    EXPECT_MATRIX_FLOAT_EQ(x, x_i);
  }
}
TEST(prob_transform, cov_matrix_constrain_exception) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, 1> x(7);
  int K = 12;
  EXPECT_THROW(stan::math::cov_matrix_constrain(x, K), std::invalid_argument);
}
TEST(prob_transform, cov_matrix_free_exception) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> y(0, 0);

  EXPECT_THROW(stan::math::cov_matrix_free(y), std::invalid_argument);
  y.resize(0, 10);
  EXPECT_THROW(stan::math::cov_matrix_free(y), std::invalid_argument);
  y.resize(10, 0);
  EXPECT_THROW(stan::math::cov_matrix_free(y), std::invalid_argument);
  y.resize(1, 2);
  EXPECT_THROW(stan::math::cov_matrix_free(y), std::invalid_argument);
  y.resize(2, 2);
  y << 0, 0, 0, 0;
  EXPECT_THROW(stan::math::cov_matrix_free(y), std::domain_error);
  y << 1, 0, 0, -1;
  EXPECT_THROW(stan::math::cov_matrix_free(y), std::domain_error);
}
TEST(covMatrixTransform, symmetry) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::cov_matrix_constrain;
  for (int K = 1; K <= 50; ++K) {
    int N = (K * (K + 1)) / 2;
    VectorXd v(N);
    for (int n = 0; n < N; ++n)
      v(n) = (n - 0.5 * N) / 10;
    MatrixXd Sigma = cov_matrix_constrain(v, K);
    EXPECT_EQ(K, Sigma.rows());
    EXPECT_EQ(K, Sigma.cols());
    for (int i = 1; i < K; ++i)
      for (int j = 0; j < i; ++j)
        EXPECT_EQ(Sigma(i, j), Sigma(j, i));  // hard equality
  }
}
