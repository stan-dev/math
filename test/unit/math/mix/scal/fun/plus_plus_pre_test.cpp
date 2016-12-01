#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

typedef stan::math::fvar<AVAR> fv_t;
typedef stan::math::fvar<fv_t> ffv_t;

TEST(MathFunctions, plusPlusPreFvarValue) {
  using stan::math::plus_plus_pre;
  fv_t x = 4;
  fv_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(5, y.val_.val());
  EXPECT_FLOAT_EQ(5, x.val_.val());
}

TEST(MathFunctions, plusPlusPreFvarDeriv) {
  using stan::math::plus_plus_pre;
  fv_t x(4, 2.3);
  fv_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(2.3, y.d_.val());
}

TEST(MathFunctions, plusPlusPreFvarRef) {
  using stan::math::plus_plus_pre;
  fv_t x(4, 2.3);
  EXPECT_EQ(&x, &plus_plus_pre(x));
}


TEST(MathFunctions, plusPlusPreF2varValue) {
  using stan::math::plus_plus_pre;
  ffv_t x(4, 2.3);
  ffv_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(5, y.val_.val_.val());
  EXPECT_FLOAT_EQ(5, x.val_.val_.val());
}

TEST(MathFunctions, plusPlusPreF2varDeriv) {
  using stan::math::plus_plus_pre;
  ffv_t x(4, 2.3);
  ffv_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(2.3, y.d_.val_.val());
  EXPECT_FLOAT_EQ(2.3, x.d_.val_.val());
}

TEST(MathFunctions, plusPlusPreF2varRef) {
  using stan::math::plus_plus_pre;
  ffv_t x(4, 2.3);
  EXPECT_EQ(&x, &plus_plus_pre(x));
}
