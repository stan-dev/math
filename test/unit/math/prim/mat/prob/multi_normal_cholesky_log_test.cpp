#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbMultiNormalCholesky, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L
      = Sigma.llt().matrixL();

  EXPECT_FLOAT_EQ((stan::math::multi_normal_cholesky_lpdf(y, mu, L)),
                  (stan::math::multi_normal_cholesky_log(y, mu, L)));
  EXPECT_FLOAT_EQ((stan::math::multi_normal_cholesky_lpdf<true>(y, mu, L)),
                  (stan::math::multi_normal_cholesky_log<true>(y, mu, L)));
  EXPECT_FLOAT_EQ((stan::math::multi_normal_cholesky_lpdf<false>(y, mu, L)),
                  (stan::math::multi_normal_cholesky_log<false>(y, mu, L)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_normal_cholesky_lpdf<true, Eigen::VectorXd,
                                              Eigen::VectorXd, Eigen::MatrixXd>(
          y, mu, L)),
      (stan::math::multi_normal_cholesky_log<true, Eigen::VectorXd,
                                             Eigen::VectorXd, Eigen::MatrixXd>(
          y, mu, L)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_normal_cholesky_lpdf<false, Eigen::VectorXd,
                                              Eigen::VectorXd, Eigen::MatrixXd>(
          y, mu, L)),
      (stan::math::multi_normal_cholesky_log<false, Eigen::VectorXd,
                                             Eigen::VectorXd, Eigen::MatrixXd>(
          y, mu, L)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_normal_cholesky_lpdf<Eigen::VectorXd, Eigen::VectorXd,
                                              Eigen::MatrixXd>(y, mu, L)),
      (stan::math::multi_normal_cholesky_log<Eigen::VectorXd, Eigen::VectorXd,
                                             Eigen::MatrixXd>(y, mu, L)));
}
