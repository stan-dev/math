#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdLogInvLogitDiff, Fvar) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;

  fvar<double> x(0.5, 1.0);
  fvar<double> y(-1.0, 1.0);
  fvar<double> z(0.0, 3.0);

  fvar<double> a = log_inv_logit_diff(x, y);
  EXPECT_FLOAT_EQ(-1.039821131, a.val_);
  EXPECT_FLOAT_EQ(0.1085992474, a.d_);

  fvar<double> b = log_inv_logit_diff(x, z);
  EXPECT_FLOAT_EQ(-2.09997629431, b.val_);
  EXPECT_FLOAT_EQ(-4.20544749628, b.d_);
}
TEST(AgradFwdLogInvLogit, FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;

  fvar<fvar<double> > x;
  x.val_.val_ = 3.4;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > y;
  y.val_.val_ = 0.9;
  y.val_.d_ = 1.0;

  fvar<fvar<double> > a = log_inv_logit_diff(x, y);

  EXPECT_FLOAT_EQ(-1.3596328289, a.val_.val_);
  EXPECT_FLOAT_EQ(-0.678654037927, a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);
}
