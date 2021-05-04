#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbMultiNormal, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  EXPECT_FLOAT_EQ((stan::math::multi_normal_lpdf(y, mu, Sigma)),
                  (stan::math::multi_normal_log(y, mu, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::multi_normal_lpdf<true>(y, mu, Sigma)),
                  (stan::math::multi_normal_log<true>(y, mu, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::multi_normal_lpdf<false>(y, mu, Sigma)),
                  (stan::math::multi_normal_log<false>(y, mu, Sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_normal_lpdf<true, Eigen::VectorXd, Eigen::VectorXd,
                                     Eigen::MatrixXd>(y, mu, Sigma)),
      (stan::math::multi_normal_log<true, Eigen::VectorXd, Eigen::VectorXd,
                                    Eigen::MatrixXd>(y, mu, Sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_normal_lpdf<false, Eigen::VectorXd, Eigen::VectorXd,
                                     Eigen::MatrixXd>(y, mu, Sigma)),
      (stan::math::multi_normal_log<false, Eigen::VectorXd, Eigen::VectorXd,
                                    Eigen::MatrixXd>(y, mu, Sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_normal_lpdf<Eigen::VectorXd, Eigen::VectorXd,
                                     Eigen::MatrixXd>(y, mu, Sigma)),
      (stan::math::multi_normal_log<Eigen::VectorXd, Eigen::VectorXd,
                                    Eigen::MatrixXd>(y, mu, Sigma)));
}
