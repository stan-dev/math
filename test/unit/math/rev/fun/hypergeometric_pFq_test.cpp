#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST(RevMath, hypergeometric_pFq) {
  using stan::math::hypergeometric_pFq;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_v a(2);
  a << 4, 2;
  vector_v b(2);
  b << 6, 3;
  var z = 4;

  var result = hypergeometric_pFq(a, b, z);
  result.grad();

  EXPECT_FLOAT_EQ(3.924636646666071, a[0].adj());
  EXPECT_FLOAT_EQ(6.897245961898751, a[1].adj());

  EXPECT_FLOAT_EQ(-2.775051002566842, b[0].adj());
  EXPECT_FLOAT_EQ(-4.980095849781222, b[1].adj());

  EXPECT_FLOAT_EQ(4.916522138006060, z.adj());
}
