#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbMultiGpCholesky, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu(5, 1);
  mu.setZero();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L
      = Sigma.llt().matrixL();

  Eigen::Matrix<double, Eigen::Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  EXPECT_FLOAT_EQ((stan::math::multi_gp_cholesky_lpdf(y, L, w)),
                  (stan::math::multi_gp_cholesky_log(y, L, w)));
  EXPECT_FLOAT_EQ((stan::math::multi_gp_cholesky_lpdf<true>(y, L, w)),
                  (stan::math::multi_gp_cholesky_log<true>(y, L, w)));
  EXPECT_FLOAT_EQ((stan::math::multi_gp_cholesky_lpdf<false>(y, L, w)),
                  (stan::math::multi_gp_cholesky_log<false>(y, L, w)));
  EXPECT_FLOAT_EQ((stan::math::multi_gp_cholesky_lpdf<true>(y, L, w)),
                  (stan::math::multi_gp_cholesky_log<true>(y, L, w)));
  EXPECT_FLOAT_EQ((stan::math::multi_gp_cholesky_lpdf<false>(y, L, w)),
                  (stan::math::multi_gp_cholesky_log<false>(y, L, w)));
  EXPECT_FLOAT_EQ((stan::math::multi_gp_cholesky_lpdf<>(y, L, w)),
                  (stan::math::multi_gp_cholesky_log<>(y, L, w)));
}
