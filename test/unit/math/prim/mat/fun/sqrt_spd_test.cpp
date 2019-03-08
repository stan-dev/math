#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>

TEST(MathMatrix, sqrt_spd) {
  stan::math::matrix_d m0;
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_d ev_m1(1, 1);
  ev_m1 << 4.0;

  using stan::math::sqrt_spd;
  EXPECT_THROW(sqrt_spd(m0), std::invalid_argument);
  EXPECT_NEAR(2.0, sqrt_spd(ev_m1)(0, 0), 1e-16);
  EXPECT_THROW(sqrt_spd(m1), std::invalid_argument);
  ev_m1(0, 0) = -1;
  EXPECT_TRUE(std::isnan(sqrt_spd(ev_m1)(0, 0)));
}
