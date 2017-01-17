#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbMultiStudentT, log_matches_lpmf) {
  Eigen::Matrix<double,Eigen::Dynamic,1> y(3,1);
  y << 2.0, -2.0, 11.0;
  Eigen::Matrix<double,Eigen::Dynamic,1> mu(3,1);
  mu << 1.0, -1.0, 3.0;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Sigma(3,3);
  Sigma << 9.0, -3.0, 0.0,
    -3.0,  4.0, 0.0,
    0.0, 0.0, 5.0;
  double nu = 4.0;
  EXPECT_FLOAT_EQ((stan::math::multi_student_t_lpdf(y,nu,mu,Sigma)),
                  (stan::math::multi_student_t_log(y,nu,mu,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::multi_student_t_lpdf<true>(y,nu,mu,Sigma)),
                  (stan::math::multi_student_t_log<true>(y,nu,mu,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::multi_student_t_lpdf<false>(y,nu,mu,Sigma)),
                  (stan::math::multi_student_t_log<false>(y,nu,mu,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::multi_student_t_lpdf<true, Eigen::VectorXd, double, Eigen::VectorXd, Eigen::MatrixXd>(y,nu,mu,Sigma)),
                  (stan::math::multi_student_t_log<true, Eigen::VectorXd, double, Eigen::VectorXd, Eigen::MatrixXd>(y,nu,mu,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::multi_student_t_lpdf<false, Eigen::VectorXd, double, Eigen::VectorXd, Eigen::MatrixXd>(y,nu,mu,Sigma)),
                  (stan::math::multi_student_t_log<false, Eigen::VectorXd, double, Eigen::VectorXd, Eigen::MatrixXd>(y,nu,mu,Sigma)));
  EXPECT_FLOAT_EQ((stan::math::multi_student_t_lpdf<Eigen::VectorXd, double, Eigen::VectorXd, Eigen::MatrixXd>(y,nu,mu,Sigma)),
                  (stan::math::multi_student_t_log<Eigen::VectorXd, double, Eigen::VectorXd, Eigen::MatrixXd>(y,nu,mu,Sigma)));
}
