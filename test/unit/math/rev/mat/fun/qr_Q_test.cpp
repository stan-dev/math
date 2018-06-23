#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(MathMatrix, qr_Q) {
  stan::math::matrix_v m0(0, 0);
  stan::math::matrix_v m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;

  using stan::math::qr_Q;
  using stan::math::transpose;
  EXPECT_THROW(qr_Q(m0), std::invalid_argument);
  EXPECT_NO_THROW(qr_Q(m1));
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::matrix_v m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;
  test::check_varis_on_stack(stan::math::qr_Q(m1));
}
