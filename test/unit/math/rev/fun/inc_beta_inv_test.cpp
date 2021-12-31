#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(inc_beta_inv, values) {
  using stan::math::var;
  using stan::math::inc_beta_inv;

  var a = 25;
  var b = 2;
  var p = 0.7;

  var res = inc_beta_inv(a, b, p);
  res.grad();

  EXPECT_FLOAT_EQ(a.adj(), 0.00161775443606);
  EXPECT_FLOAT_EQ(b.adj(), -0.0289198999936);
  EXPECT_FLOAT_EQ(p.adj(), 0.102597843424);
}
