#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, qr_thin_Q) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;

  using stan::math::qr_thin_Q;
  using stan::math::transpose;
  EXPECT_NO_THROW(qr_thin_Q(m0));
  EXPECT_NO_THROW(qr_thin_Q(m1));
  EXPECT_NO_THROW(qr_thin_Q(transpose(m1)));
}
