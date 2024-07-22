#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, qr_thin) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(4, 2);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8;

  using stan::math::qr_thin;
  using stan::math::transpose;
  EXPECT_NO_THROW(qr_thin(m0));
  EXPECT_NO_THROW(qr_thin(m1));

  stan::math::matrix_d m2(4, 2);
  stan::math::matrix_d Q(4, 2);
  stan::math::matrix_d R(2, 2);
  std::tie(Q, R) = qr_thin(m1);
  m2 = Q * R;
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      EXPECT_NEAR(m1(i, j), m2(i, j), 1e-12);
    }
  }
  for (unsigned int j = 0; j < m1.cols(); j++)
    EXPECT_TRUE(R(j, j) >= 0.0);

  auto m1_T = transpose(m1);

  std::tie(Q, R) = qr_thin(m1_T);
  stan::math::matrix_d m3(2, 4);
  m3 = Q * R;
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      EXPECT_NEAR(m1(i, j), m3(j, i), 1e-12);
    }
  }
}
