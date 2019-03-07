#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>

TEST(MathMatrix, inv_sqrt_spd) {
  stan::math::matrix_d m0;
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_d ev_m1(1, 1);
  ev_m1 << 4.0;

  using stan::math::inv_sqrt_spd;
  EXPECT_THROW(inv_sqrt_spd(m0), std::invalid_argument);
  EXPECT_NEAR(0.5, inv_sqrt_spd(ev_m1)(0, 0), 1e-16);
  EXPECT_THROW(inv_sqrt_spd(m1), std::invalid_argument);
}
