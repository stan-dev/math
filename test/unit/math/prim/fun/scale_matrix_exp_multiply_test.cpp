#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <vector>

template <int N, int M>
inline void test_scale_matrix_exp_multiply() {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::scale_matrix_exp_multiply;

  std::srand(1999);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(N, N);
  Eigen::Matrix<double, -1, M> B = Eigen::MatrixXd::Random(N, M);
  Eigen::Matrix<double, Dynamic, Dynamic> A0 = A;

  // brute force
  Eigen::Matrix<double, N, M> expAB
      = stan::math::multiply(stan::math::matrix_exp(A0), B);

  // matrix_exp_multiply
  const double t = 1.0;
  Eigen::Matrix<double, N, M> res = scale_matrix_exp_multiply(t, A, B);
  EXPECT_MATRIX_FLOAT_EQ(expAB, res);
}

TEST(MathMatrixPrimMat, scale_matrix_exp_multiply) {
  // the helper above doesn't handle 0 size inputs
  const double t = 1.0;
  Eigen::MatrixXd A(0, 0);
  Eigen::MatrixXd B(0, 0);
  EXPECT_EQ(stan::math::scale_matrix_exp_multiply(t, A, B).size(), 0);

  Eigen::MatrixXd C(0, 2);
  Eigen::MatrixXd M = stan::math::scale_matrix_exp_multiply(t, A, C);
  EXPECT_EQ(A.rows(), M.rows());
  EXPECT_EQ(C.cols(), M.cols());

  test_scale_matrix_exp_multiply<1, 1>();
  test_scale_matrix_exp_multiply<1, 5>();
  test_scale_matrix_exp_multiply<5, 1>();
  test_scale_matrix_exp_multiply<5, 5>();
  test_scale_matrix_exp_multiply<20, 2>();
}

TEST(MathMatrixPrimMat, scale_matrix_exp_multiply_issue_1146) {
  // matrix from bug report
  // https://github.com/stan-dev/math/issues/1146
  int m = 17;
  Eigen::MatrixXd Ac(m, m);
  Ac << -4.4758400, 0.0732692, 0.0000000, 0.000000, 0.000000, 0.0000000, 0.0000,
      0.00000000, 0.000000, 0.0000000, 0.00000000, 0.0000000, 0.0000000,
      0.0000000, 0.00000000, 0.00000000, 0.00000000, 0.0183173, -0.0366346,
      0.0183173, 0.000000, 0.000000, 0.0000000, 0.0000, 0.00000000, 0.000000,
      0.0000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.00000000,
      0.00000000, 0.00000000, 0.0000000, 0.0146538, -0.1794070, 0.000000,
      0.000000, 0.0268619, 0.0000, 0.00639568, 0.000000, 0.0639568, 0.00000000,
      0.0000000, 0.0426379, 0.0249012, 0.00748033, 0.00498689, 0.00415574,
      0.0000000, 0.0000000, 0.0000000, -5.868960, 0.447212, 0.0000000, 0.0000,
      0.00000000, 0.000000, 0.0000000, 0.00000000, 0.0000000, 0.0000000,
      0.0000000, 0.00000000, 0.00000000, 0.00000000, 0.0000000, 0.0000000,
      0.0000000, 0.223606, -0.447212, 0.2236060, 0.0000, 0.00000000, 0.000000,
      0.0000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.00000000,
      0.00000000, 0.00000000, 0.0000000, 0.0000000, 0.3287580, 0.000000,
      0.447212, -2.3328900, 0.0000, 0.04931380, 0.000000, 0.4931380, 0.00000000,
      0.0000000, 0.3287580, 0.6857140, 0.05767690, 0.03845130, 0.03204270,
      0.0000000, 0.0000000, 0.0000000, 0.000000, 0.000000, 0.0000000, -42.1500,
      11.43980000, 0.000000, 0.0000000, 0.00000000, 0.0000000, 0.0000000,
      0.0000000, 0.00000000, 0.00000000, 0.00000000, 0.0000000, 0.0000000,
      2.3972000, 0.000000, 0.000000, 1.5102300, 11.4398, -26.34030000, 0.000000,
      3.5957900, 0.00000000, 0.0000000, 2.3972000, 5.0000000, 0.42056100,
      0.28037400, 0.23364500, 0.0000000, 0.0000000, 0.0000000, 0.000000,
      0.000000, 0.0000000, 0.0000, 0.00000000, -0.122394, 0.1223940, 0.00000000,
      0.0000000, 0.0000000, 0.0000000, 0.00000000, 0.00000000, 0.00000000,
      0.0000000, 0.0000000, 1.0958600, 0.000000, 0.000000, 0.6903930, 0.0000,
      0.16437900, 1.101550, -6.4337600, 0.00000000, 0.0000000, 1.0958600,
      2.2857100, 0.19225600, 0.12817100, 0.10680900, 0.0000000, 0.0000000,
      0.0000000, 0.000000, 0.000000, 0.0000000, 0.0000, 0.00000000, 0.000000,
      0.0000000, -0.00471794, 0.0014969, 0.0000000, 0.0000000, 0.00000000,
      0.00000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.000000,
      0.000000, 0.0000000, 0.0000, 0.00000000, 0.000000, 0.0000000, 0.01901160,
      -0.0402984, 0.0212869, 0.0000000, 0.00000000, 0.00000000, 0.00000000,
      0.0000000, 0.0000000, 0.0426379, 0.000000, 0.000000, 0.0268619, 0.0000,
      0.00639568, 0.000000, 0.0639568, 0.00000000, 0.0170295, -0.3347470,
      0.1778660, 0.00748033, 0.00498689, 0.00415574, 0.0000000, 0.0000000,
      0.2695880, 0.000000, 0.000000, 0.6065730, 0.0000, 0.14442200, 0.000000,
      1.4442200, 0.00000000, 0.0000000, 1.9256300, -4.3904300, 0.03851260,
      0.15405000, 0.19256300, 0.0000000, 0.0000000, 0.0000000, 0.000000,
      0.000000, 0.0000000, 0.0000, 0.00000000, 0.000000, 0.0000000, 0.00000000,
      0.0000000, 0.0000000, 0.0000000, -1.00000000, 0.00000000, 0.00000000,
      0.0000000, 0.0000000, 0.0000000, 0.000000, 0.000000, 0.0000000, 0.0000,
      0.00000000, 0.000000, 0.0000000, 0.00000000, 0.0000000, 0.0000000,
      0.0000000, 0.00000000, -1.00000000, 0.00000000, 0.0000000, 0.0000000,
      0.0000000, 0.000000, 0.000000, 0.0000000, 0.0000, 0.00000000, 0.000000,
      0.0000000, 0.00000000, 0.0000000, 0.0000000, 0.0000000, 0.00000000,
      0.00000000, -1.00000000;

  Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(2 * m, 2 * m);
  Eigen::MatrixXd OI = Eigen::MatrixXd::Zero(2 * m, m);

  OI.bottomRightCorner(m, m) = Eigen::VectorXd::Constant(m, 1.0).asDiagonal();
  phi.topLeftCorner(m, m) = Ac;
  phi.topRightCorner(m, m) = Eigen::VectorXd::Constant(m, 0.1).asDiagonal();
  phi.bottomRightCorner(m, m) = -Ac.transpose();

  Eigen::MatrixXd expAB(2 * m, m);
  Eigen::MatrixXd res(2 * m, m);
  for (auto i = 0; i < 20; ++i) {
    double dt = 0.1 + i * 0.05;
    expAB = stan::math::matrix_exp(phi * dt) * OI;
    res = stan::math::scale_matrix_exp_multiply(dt, phi, OI);
    EXPECT_MATRIX_FLOAT_EQ(expAB, res);
  }
}

TEST(MathMatrixPrimMat, scale_matrix_exp_multiply_exception) {
  using stan::math::scale_matrix_exp_multiply;
  const double t = 1.0;
  {  // multiplicable
    Eigen::MatrixXd A(0, 0);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(1, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
    EXPECT_THROW(scale_matrix_exp_multiply(t, B, A), std::invalid_argument);
  }

  {  // multiplicable
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 2);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
  }

  {  // square
    Eigen::MatrixXd A(0, 1);
    Eigen::MatrixXd B(1, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
  }

  {  // square
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
  }
}
