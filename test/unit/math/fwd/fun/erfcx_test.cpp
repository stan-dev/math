#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(AgradFwdErfcx, Fvar) {
  using stan::math::erfcx;
  using stan::math::fvar;

  for (double x_value : {-2.0, -1.0, 0.0, 0.5, 2.0, 12.0}) {
    fvar<double> x(x_value, 1.0);
    fvar<double> y = erfcx(x);
    EXPECT_FLOAT_EQ(erfcx(x_value), y.val_);
    EXPECT_FLOAT_EQ(
        2.0 * x_value * erfcx(x_value) - stan::math::TWO_OVER_SQRT_PI, y.d_);
  }
}

TEST(AgradFwdErfcx, FvarChainRule) {
  using stan::math::erfcx;
  using stan::math::fvar;
  fvar<double> x(1.5, 3.0);
  fvar<double> y = erfcx(x);
  EXPECT_FLOAT_EQ(3.0 * (2.0 * 1.5 * erfcx(1.5) - stan::math::TWO_OVER_SQRT_PI),
                  y.d_);
}

TEST(AgradFwdErfcx, FvarFvarDouble) {
  using stan::math::erfcx;
  using stan::math::fvar;

  const double x_value = 0.75;
  fvar<fvar<double> > x;
  x.val_.val_ = x_value;
  x.val_.d_ = 1.0;
  x.d_.val_ = 1.0;

  fvar<fvar<double> > y = erfcx(x);
  const double v = erfcx(x_value);
  const double d1 = 2.0 * x_value * v - stan::math::TWO_OVER_SQRT_PI;
  // d2 = 2 * erfcx + 2 * x * d1
  const double d2 = 2.0 * v + 2.0 * x_value * d1;

  EXPECT_FLOAT_EQ(v, y.val_.val_);
  EXPECT_FLOAT_EQ(d1, y.val_.d_);
  EXPECT_FLOAT_EQ(d1, y.d_.val_);
  EXPECT_FLOAT_EQ(d2, y.d_.d_);
}
