#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, prod) {
  using stan::math::prod;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d vd;
  vector_v vv;
  EXPECT_FLOAT_EQ(1.0, prod(vd));
  EXPECT_FLOAT_EQ(1.0, prod(vv).val());

  vd = vector_d(1);
  vv = vector_v(1);
  vd << 2.0;
  vv << 2.0;
  EXPECT_FLOAT_EQ(2.0, prod(vd));
  EXPECT_FLOAT_EQ(2.0, prod(vv).val());

  vd = vector_d(2);
  vd << 2.0, 3.0;
  vv = vector_v(2);
  vv << 2.0, 3.0;
  AVEC x(2);
  x[0] = vv[0];
  x[1] = vv[1];
  AVAR f = prod(vv);
  EXPECT_FLOAT_EQ(6.0, prod(vd));
  EXPECT_FLOAT_EQ(6.0, f.val());
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(3.0, g[0]);
  EXPECT_FLOAT_EQ(2.0, g[1]);
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::vector_v v(1);
  v << 2.0;
  test::check_varis_on_stack(stan::math::prod(v));
}
