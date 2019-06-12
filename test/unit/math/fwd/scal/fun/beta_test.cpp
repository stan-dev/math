#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdBeta, Fvar) {
  using stan::math::beta;
  using stan::math::fvar;
  using stan::math::lbeta;
  using std::exp;

  fvar<double> x(0.5, 2.5);
  fvar<double> y(3.3, 3.0);
  double x_dbl = 0.5;
  double y_dbl = 3.3;

  fvar<double> a = beta(x, y);
  EXPECT_FLOAT_EQ(1.0132280381701558382291101, a.val_);
  EXPECT_FLOAT_EQ(-8.502606179100956692649726, a.d_);

  fvar<double> b = exp(lbeta(x, y));
  EXPECT_FLOAT_EQ(a.val_, b.val_);
  EXPECT_FLOAT_EQ(a.d_, b.d_);

  a = beta(x, y_dbl);
  EXPECT_FLOAT_EQ(1.0132280381701558382291101, a.val_);
  EXPECT_FLOAT_EQ(-8.0075408614604666560873932, a.d_);

  a = beta(x_dbl, y);
  EXPECT_FLOAT_EQ(1.0132280381701558382291101, a.val_);
  EXPECT_FLOAT_EQ(-0.4950653176404900365623333, a.d_);
}

TEST(AgradFwdBeta, FvarFvarDouble) {
  using stan::math::beta;
  using stan::math::fvar;
  using stan::math::lbeta;
  using std::exp;

  fvar<fvar<double> > x;
  x.val_.val_ = 3.4;
  x.val_.d_ = 5.1;

  fvar<fvar<double> > y;
  y.val_.val_ = 0.9;
  y.val_.d_ = 3.7;

  double x_dbl = 3.4;
  double y_dbl = 0.9;

  fvar<fvar<double> > a = beta(x, y);
  EXPECT_FLOAT_EQ(0.35976049995522196968654, a.val_.val_);
  EXPECT_FLOAT_EQ(-3.27797158126448572143610, a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);

  fvar<fvar<double> > b = exp(lbeta(x, y));
  EXPECT_FLOAT_EQ(a.val_.val_, b.val_.val_);
  EXPECT_FLOAT_EQ(a.val_.d_, b.val_.d_);
  EXPECT_FLOAT_EQ(a.d_.val_, b.d_.val_);
  EXPECT_FLOAT_EQ(a.d_.d_, b.d_.d_);

  a = beta(x, y_dbl);
  EXPECT_FLOAT_EQ(0.35976049995522196968654, a.val_.val_);
  EXPECT_FLOAT_EQ(-0.49224348198636858027266, a.val_.d_);

  a = beta(x_dbl, y);
  EXPECT_FLOAT_EQ(0.35976049995522196968654, a.val_.val_);
  EXPECT_FLOAT_EQ(-2.78572809927811714116343, a.val_.d_);
}
