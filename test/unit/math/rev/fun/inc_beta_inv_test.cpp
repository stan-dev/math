#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(inc_beta_inv, values) {
  using stan::math::var;
  using stan::math::inc_beta_inv;

  // check autodiffing works with var types with large values
  var a = 1;
  var b = 2;
  var p = 0.5;

  var res = inc_beta_inv(a, b, p);
  res.grad();

  EXPECT_FLOAT_EQ(a.adj(), 0.287698278597);
  EXPECT_FLOAT_EQ(b.adj(), -0.122532267934);
  EXPECT_FLOAT_EQ(p.adj(), 0.707106781187);
}
