#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, qr_R) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(4, 2);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8;

  using stan::math::qr_Q;
  using stan::math::qr_R;
  using stan::math::transpose;
  EXPECT_THROW(qr_R(m0), std::invalid_argument);
  EXPECT_NO_THROW(qr_R(m1));

  stan::math::matrix_d m2(4, 2);
  stan::math::matrix_d Q(4, 4);
  stan::math::matrix_d R(4, 2);
  Q = qr_Q(m1);
  R = qr_R(m1);
  m2 = Q * R;
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      EXPECT_NEAR(m1(i, j), m2(i, j), 1e-12);
    }
  }
  for (unsigned int j = 0; j < m1.cols(); j++)
    EXPECT_TRUE(R(j, j) >= 0.0);
  stan::math::matrix_d m3(2, 4);
  m3 = qr_Q(transpose(m1)) * qr_R(transpose(m1));
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      EXPECT_NEAR(m1(i, j), m3(j, i), 1e-12);
    }
  }
}
