#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbMatrixNormalPrec, log_matches_lpmf) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mu(3,5);
  mu.setZero();
  
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y(3,5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0,
       11.0, 2.0, -5.0, 11.0, 0.0,
       -2.0, 11.0, 2.0, -2.0, -11.0;

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Sigma(5,5);
  Sigma << 9.0, -3.0, 0.0,  0.0, 0.0,
          -3.0,  4.0, 0.0,  0.0, 0.0,
           0.0,  0.0, 5.0,  1.0, 0.0,
           0.0,  0.0, 1.0, 10.0, 0.0,
           0.0,  0.0, 0.0,  0.0, 2.0;

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> D(3,3);
  D << 1.0, 0.5, 0.1,
       0.5, 1.0, 0.2,
       0.1, 0.2, 1.0;
  
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf(y,mu,D,Sigma)),
                  (stan::math::matrix_normal_prec_log(y,mu,D,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf<true>(y,mu,D,Sigma)),
                  (stan::math::matrix_normal_prec_log<true>(y,mu,D,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf<false>(y,mu,D,Sigma)),
                  (stan::math::matrix_normal_prec_log<false>(y,mu,D,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf<true, double, double, double>(y,mu,D,Sigma)),
                  (stan::math::matrix_normal_prec_log<true, double, double, double>(y,mu,D,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf<false, double, double, double>(y,mu,D,Sigma)),
                  (stan::math::matrix_normal_prec_log<false, double, double, double>(y,mu,D,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::matrix_normal_prec_lpdf<double, double, double>(y,mu,D,Sigma)),
                  (stan::math::matrix_normal_prec_log<double, double, double>(y,mu,D,Sigma)));

}
