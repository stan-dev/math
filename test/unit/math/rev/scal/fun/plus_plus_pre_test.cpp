#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(MathFunctions, plusPlusPreVarValue) {
  using stan::math::plus_plus_pre;
  AVAR x = 4;
  AVAR y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(5, y.val());
  EXPECT_FLOAT_EQ(5, x.val());
}

TEST(MathFunctions, plusPlusPreVarDeriv) {
  using stan::math::plus_plus_pre;
  AVAR x = 4;
  AVAR y = plus_plus_pre(x);
  AVEC xs = createAVEC(x);
  VEC g;
  y.grad(xs, g);
  EXPECT_FLOAT_EQ(1, g[0]);
}

TEST(MathFunctions, plusPlusPreVarRef) {
  using stan::math::plus_plus_pre;
  AVAR x = 4;
  EXPECT_EQ(&x, &plus_plus_pre(x));
}
