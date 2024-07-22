#ifndef TEST_UNIT_MATH_PRIM_PROB_AGRAD_DISTRIBUTIONS_MULTI_STUDENT_T_CHOLESKY
#define TEST_UNIT_MATH_PRIM_PROB_AGRAD_DISTRIBUTIONS_MULTI_STUDENT_T_CHOLESKY
class agrad_distributions_multi_student_t_cholesky : public ::testing::Test {
 protected:
  virtual void SetUp() {
    nu = 5;

    y.resize(3, 1);
    y << 2.0, -2.0, 11.0;
    y2.resize(3, 1);
    y2 << 15.0, 1.0, -5.0;

    mu.resize(3, 1);
    mu << 1.0, -1.0, 3.0;
    mu2.resize(3, 1);
    mu2 << 6.0, 2.0, -6.0;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma(3, 3);
    Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
    L = Sigma.llt().matrixL();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma2(3, 3);
    Sigma2 << 3.0, 1.0, 0.0, 1.0, 5.0, -2.0, 0.0, -2.0, 9.0;
    L2 = Sigma2.llt().matrixL();
  }

  double nu;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu;
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu2;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L2;
};
#endif
