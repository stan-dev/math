#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdBeta, Fvar) {
  using stan::math::fvar;
  using stan::math::beta;
  using stan::math::lbeta;
  using std::exp;

  fvar<double> x(0.5, 1.0);
  fvar<double> y(3.3, 1.0);

  fvar<double> a = beta(x, y);
  EXPECT_FLOAT_EQ(1.0132280381701558382291101, a.val_);
  EXPECT_FLOAT_EQ(-3.368038117131016674622401, a.d_);

  fvar<double> b = exp(lbeta(x, y));
  EXPECT_FLOAT_EQ(a.val_, b.val_);
  EXPECT_FLOAT_EQ(a.d_, b.d_);
}

TEST(AgradFwdBeta, FvarFvarDouble) {
  using stan::math::fvar;
  using stan::math::beta;
  using stan::math::lbeta;
  using std::exp;

  fvar<fvar<double> > x;
  x.val_.val_ = 3.4;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > y;
  y.val_.val_ = 0.9;
  y.val_.d_ = 1.0;

  fvar<fvar<double> > a = beta(x, y);
  EXPECT_FLOAT_EQ(0.35976049995522196968654, a.val_.val_);
  EXPECT_FLOAT_EQ(-0.8494178160926317523551, a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);

  fvar<fvar<double> > b = exp(lbeta(x, y));
  EXPECT_FLOAT_EQ(a.val_.val_, b.val_.val_);
  EXPECT_FLOAT_EQ(a.val_.d_, b.val_.d_);
  EXPECT_FLOAT_EQ(a.d_.val_, b.d_.val_);
  EXPECT_FLOAT_EQ(a.d_.d_, b.d_.d_);
}
