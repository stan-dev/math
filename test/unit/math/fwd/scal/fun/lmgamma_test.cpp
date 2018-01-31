#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdLmgamma, Fvar) {
  using stan::math::fvar;
  using stan::math::lmgamma;

  int x = 3;
  fvar<double> y(3.2, 2.1);

  fvar<double> a = lmgamma(x, y);
  EXPECT_FLOAT_EQ(lmgamma(3, 3.2), a.val_);
  EXPECT_FLOAT_EQ(4.9138227, a.d_);
}

TEST(AgradFwdLmgamma, FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::lmgamma;

  fvar<fvar<double> > x;
  x.val_.val_ = 3.2;
  x.val_.d_ = 2.1;

  fvar<fvar<double> > a = lmgamma(3, x);

  EXPECT_FLOAT_EQ(lmgamma(3, 3.2), a.val_.val_);
  EXPECT_FLOAT_EQ(4.9138227, a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);

  fvar<fvar<double> > y;
  y.val_.val_ = 3.2;
  y.d_.val_ = 2.1;

  a = lmgamma(3, y);
  EXPECT_FLOAT_EQ(lmgamma(3, 3.2), a.val_.val_);
  EXPECT_FLOAT_EQ(0, a.val_.d_);
  EXPECT_FLOAT_EQ(4.9138227, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);
}

struct lmgamma_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return lmgamma(3, arg1);
  }
};

TEST(AgradFwdLmgamma, lmgamma_NaN) {
  lmgamma_fun lmgamma_;
  test_nan_fwd(lmgamma_, false);
}
