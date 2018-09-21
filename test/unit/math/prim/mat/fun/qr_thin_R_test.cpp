#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, qr_thin_R) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(4, 2);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8;

  using stan::math::qr_thin_Q;
  using stan::math::qr_thin_R;
  using stan::math::transpose;
  EXPECT_THROW(qr_thin_R(m0), std::invalid_argument);
  EXPECT_NO_THROW(qr_thin_R(m1));

  stan::math::matrix_d m2(4, 2);
  stan::math::matrix_d Q(4, 2);
  stan::math::matrix_d R(2, 2);
  Q = qr_thin_Q(m1);
  R = qr_thin_R(m1);
  m2 = Q * R;
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      EXPECT_NEAR(m1(i, j), m2(i, j), 1e-12);
    }
  }
  for (unsigned int j = 0; j < m1.cols(); j++)
    EXPECT_TRUE(R(j, j) >= 0.0);

  stan::math::matrix_d m3(2, 4);
  m3 = qr_thin_Q(transpose(m1)) * qr_thin_R(transpose(m1));
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      EXPECT_NEAR(m1(i, j), m3(j, i), 1e-12);
    }
  }

  stan::math::matrix_d m4(3, 3);
  m4 << -1, -2, -3, -4, -5, -6, -7, -8, -9;
  stan::math::matrix_d m5(3, 3);
  m5 = qr_thin_Q(m4) * qr_thin_R(m4);
  for (unsigned int i = 0; i < m4.rows(); i++) {
    for (unsigned int j = 0; j < m4.cols(); j++) {
      EXPECT_NEAR(m4(i, j), m5(i, j), 1e-12);
    }
  }
}
