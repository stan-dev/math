#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbMatrixNormalPrec, log_matches_lpmf) {
  using Eigen::MatrixXd;

  MatrixXd y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  MatrixXd mu = MatrixXd::Zero(3, 5);

  MatrixXd Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  MatrixXd D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D)),
                  (stan::math::matrix_normal_prec_log(y, mu, Sigma, D)));
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf<true>(y, mu, Sigma, D)),
                  (stan::math::matrix_normal_prec_log<true>(y, mu, Sigma, D)));
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf<false>(y, mu, Sigma, D)),
                  (stan::math::matrix_normal_prec_log<false>(y, mu, Sigma, D)));
}
