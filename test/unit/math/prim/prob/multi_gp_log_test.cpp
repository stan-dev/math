#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbMultiGp, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu(5, 1);
  mu.setZero();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Eigen::Matrix<double, Eigen::Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  EXPECT_FLOAT_EQ((stan::math::multi_gp_lpdf(y, Sigma, w)),
                  (stan::math::multi_gp_log(y, Sigma, w)));
  EXPECT_FLOAT_EQ((stan::math::multi_gp_lpdf<true>(y, Sigma, w)),
                  (stan::math::multi_gp_log<true>(y, Sigma, w)));
  EXPECT_FLOAT_EQ((stan::math::multi_gp_lpdf<false>(y, Sigma, w)),
                  (stan::math::multi_gp_log<false>(y, Sigma, w)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_gp_lpdf<true, double, double, double>(y, Sigma, w)),
      (stan::math::multi_gp_log<true, double, double, double>(y, Sigma, w)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_gp_lpdf<false, double, double, double>(y, Sigma, w)),
      (stan::math::multi_gp_log<false, double, double, double>(y, Sigma, w)));
  EXPECT_FLOAT_EQ(
      (stan::math::multi_gp_lpdf<double, double, double>(y, Sigma, w)),
      (stan::math::multi_gp_log<double, double, double>(y, Sigma, w)));
}
