#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(AgradRev, int_step) {
  using stan::math::int_step;

  stan::math::var a(5.0);
  stan::math::var b(0.0);
  stan::math::var c(-1.0);

  EXPECT_EQ(1U, int_step(a));
  EXPECT_EQ(0U, int_step(b));
  EXPECT_EQ(0U, int_step(c));
}
