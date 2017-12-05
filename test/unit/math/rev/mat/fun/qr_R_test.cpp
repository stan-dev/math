#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(MathMatrix, qr_R) {
  stan::math::matrix_v m0(0, 0);
  stan::math::matrix_v m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;

  using stan::math::qr_R;
  using stan::math::qr_Q;
  using stan::math::transpose;
  EXPECT_THROW(qr_R(m0), std::invalid_argument);
  EXPECT_NO_THROW(qr_R(m1));

  stan::math::matrix_v m2(3, 2);
  m2 = qr_Q(m1) * qr_R(m1);
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      EXPECT_NEAR(m1(i, j).val(), m2(i, j).val(), 1e-8);
    }
  }
  EXPECT_THROW(qr_R(transpose(m1)), std::domain_error);
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::matrix_v m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;
  test::check_varis_on_stack(stan::math::qr_R(m1));
}
