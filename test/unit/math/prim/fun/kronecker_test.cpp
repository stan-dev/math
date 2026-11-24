#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <complex>

TEST(MathMatrixPrim, kronecker) {
  Eigen::MatrixXd a{{1, -2, 3}, {0.5, 4, -1}};
  Eigen::MatrixXd b{{2, 1}, {0, -3}, {5, 0.25}};
  // numpy.kron(a, b)
  Eigen::MatrixXd expected{
      {2, 1, -4, -2, 6, 3},           {0, -3, 0, 6, 0, -9},
      {5, 0.25, -10, -0.5, 15, 0.75}, {1, 0.5, 8, 4, -2, -1},
      {0, -1.5, 0, -12, 0, 3},        {2.5, 0.125, 20, 1, -5, -0.25}};
  EXPECT_MATRIX_FLOAT_EQ(expected, stan::math::kronecker(a, b));
}

TEST(MathMatrixPrim, kronecker_complex) {
  using c = std::complex<double>;
  Eigen::MatrixXcd a{{c(1, 2), c(-1, 0)}, {c(0, 0.5), c(3, 0)}};
  Eigen::MatrixXcd b{{c(2, -1)}, {c(1, 1)}};
  // numpy.kron(a, b)
  Eigen::MatrixXcd expected{{c(4, 3), c(-2, 1)},
                            {c(-1, 3), c(-1, -1)},
                            {c(0.5, 1), c(6, -3)},
                            {c(-0.5, 0.5), c(3, 3)}};
  Eigen::MatrixXcd res = stan::math::kronecker(a, b);
  ASSERT_EQ(expected.rows(), res.rows());
  ASSERT_EQ(expected.cols(), res.cols());
  for (int i = 0; i < expected.size(); ++i) {
    EXPECT_FLOAT_EQ(expected(i).real(), res(i).real());
    EXPECT_FLOAT_EQ(expected(i).imag(), res(i).imag());
  }
}

TEST(MathMatrixPrim, kronecker_empty) {
  Eigen::MatrixXd a{{1, 2}, {3, 4}};
  Eigen::MatrixXd b(0, 3);
  Eigen::MatrixXd res = stan::math::kronecker(a, b);
  EXPECT_EQ(0, res.rows());
  EXPECT_EQ(6, res.cols());
  res = stan::math::kronecker(b, a);
  EXPECT_EQ(0, res.rows());
  EXPECT_EQ(6, res.cols());
}
