#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, diag_matrix_exception) {
  stan::math::vector_d v0;

  using stan::math::diag_matrix;
  EXPECT_NO_THROW(diag_matrix(v0));
}

TEST(MathMatrixPrimMat, diag_matrix_row_vector) {
  using stan::math::diag_matrix;
  stan::math::row_vector_d rv(3);
  rv << 1, 2, 3;
  stan::math::vector_d v = rv.transpose();
  EXPECT_MATRIX_EQ(diag_matrix(v), diag_matrix(rv));
  EXPECT_NO_THROW(diag_matrix(stan::math::row_vector_d()));

  Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> crv(2);
  crv << std::complex<double>(1, 2), std::complex<double>(3, 4);
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cm
      = diag_matrix(crv);
  EXPECT_EQ(cm(0, 0), crv(0));
  EXPECT_EQ(cm(1, 1), crv(1));
  EXPECT_EQ(cm(0, 1), std::complex<double>(0, 0));
}
