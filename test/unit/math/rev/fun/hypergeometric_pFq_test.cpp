#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <limits>

TEST(primScalFun, grad_2F2) {
  using stan::math::vector_v;
  using stan::math::vector_d;
  using stan::math::var;
  using stan::math::hypergeometric_pFq;

  vector_v p(2);
  p << 4, 2;
  vector_v q(2);
  q << 6, 3;
  var z = 4; 

  var result = hypergeometric_pFq(p, q, z);
  result.grad();

  EXPECT_FLOAT_EQ(3.924636646666071, p[0].adj());
  EXPECT_FLOAT_EQ(6.897245961898751, p[1].adj());
  
  EXPECT_FLOAT_EQ(-2.775051002566842, q[0].adj());
  EXPECT_FLOAT_EQ(-4.980095849781222, q[1].adj());

  EXPECT_FLOAT_EQ(4.916522138006060, z.adj());
}
