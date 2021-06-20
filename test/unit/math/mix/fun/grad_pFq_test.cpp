#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(ProbInternalMath, grad_2F2_vvv) {
  using stan::math::grad_pFq;
  using stan::math::vector_v;
  using stan::math::var;

  vector_v p(2);
  p << 4, 2;
  vector_v q(2);
  q << 6, 3;
  var z = 4;

  vector_v grad_p(2);
  vector_v grad_q(2);
  var grad_z;

  grad_pFq(grad_p, grad_q, grad_z, p, q, z);

  EXPECT_FLOAT_EQ(3.924636646666071, grad_p[0].val());
  EXPECT_FLOAT_EQ(6.897245961898751, grad_p[1].val());

  EXPECT_FLOAT_EQ(-2.775051002566842, grad_q[0].val());
  EXPECT_FLOAT_EQ(-4.980095849781222, grad_q[1].val());

  EXPECT_FLOAT_EQ(4.916522138006060, grad_z.val());
}
