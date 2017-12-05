#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

TEST(AgradFwdLog1pExp, Fvar) {
  using stan::math::fvar;
  using stan::math::log1p_exp;
  using std::exp;

  fvar<double> x(0.5, 1.0);
  fvar<double> y(1.0, 2.0);
  fvar<double> z(2.0, 3.0);

  fvar<double> a = log1p_exp(x);
  EXPECT_FLOAT_EQ(log1p_exp(0.5), a.val_);
  EXPECT_FLOAT_EQ(exp(0.5) / (1 + exp(0.5)), a.d_);

  fvar<double> b = log1p_exp(y);
  EXPECT_FLOAT_EQ(log1p_exp(1.0), b.val_);
  EXPECT_FLOAT_EQ(2.0 * exp(1.0) / (1 + exp(1.0)), b.d_);

  fvar<double> a2 = log(1 + exp(x));
  EXPECT_FLOAT_EQ(a.d_, a2.d_);

  fvar<double> b2 = log(1 + exp(y));
  EXPECT_FLOAT_EQ(b.d_, b2.d_);
}

TEST(AgradFwdLog1pExp, FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::log1p_exp;
  using std::exp;

  fvar<fvar<double> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > a = log1p_exp(x);

  EXPECT_FLOAT_EQ(log1p_exp(0.5), a.val_.val_);
  EXPECT_FLOAT_EQ(exp(0.5) / (1 + exp(0.5)), a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);

  fvar<fvar<double> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  a = log1p_exp(y);
  EXPECT_FLOAT_EQ(log1p_exp(0.5), a.val_.val_);
  EXPECT_FLOAT_EQ(0, a.val_.d_);
  EXPECT_FLOAT_EQ(exp(0.5) / (1 + exp(0.5)), a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);
}

struct log1p_exp_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return log1p_exp(arg1);
  }
};

TEST(AgradFwdLog1pExp, log1p_exp_NaN) {
  log1p_exp_fun log1p_exp_;
  test_nan_fwd(log1p_exp_, false);
}
