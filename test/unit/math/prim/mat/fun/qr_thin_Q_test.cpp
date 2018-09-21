#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, qr_thin_Q) {
  stan::math::matrix_d m0(0, 0);
  stan::math::matrix_d m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;

  using stan::math::qr_thin_Q;
  using stan::math::transpose;
  EXPECT_THROW(qr_thin_Q(m0), std::invalid_argument);
  EXPECT_NO_THROW(qr_thin_Q(m1));
  EXPECT_NO_THROW(qr_thin_Q(transpose(m1)));
}
