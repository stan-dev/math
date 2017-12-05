#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>
#include <limits>

TEST(AgradFwdIsNan, Fvar) {
  using stan::math::fvar;
  using stan::math::is_nan;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  fvar<double> a(nan, nan);
  fvar<double> b(max, max);
  fvar<double> c(min, min);
  fvar<double> d(0.5, 1.0);
  fvar<double> e(infinity, infinity);

  EXPECT_TRUE(is_nan(a.val_));
  EXPECT_FALSE(is_nan(b.val_));
  EXPECT_FALSE(is_nan(c.val_));
  EXPECT_FALSE(is_nan(d.val_));
  EXPECT_FALSE(is_nan(e.val_));
}

TEST(AgradFwdIsNan, FvarFvar) {
  using stan::math::fvar;
  using stan::math::is_nan;

  double infinity = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  double min = std::numeric_limits<double>::min();
  double max = std::numeric_limits<double>::max();

  fvar<fvar<double> > a(nan, nan);
  fvar<fvar<double> > b(max, max);
  fvar<fvar<double> > c(min, min);
  fvar<fvar<double> > d(0.5, 1.0);
  fvar<fvar<double> > e(infinity, infinity);

  EXPECT_TRUE(is_nan(a.val_.val_));
  EXPECT_FALSE(is_nan(b.val_.val_));
  EXPECT_FALSE(is_nan(c.val_.val_));
  EXPECT_FALSE(is_nan(d.val_.val_));
  EXPECT_FALSE(is_nan(e.val_.val_));
}

