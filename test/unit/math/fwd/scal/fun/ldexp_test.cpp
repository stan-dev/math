#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/fwd/scal/fun/nan_util.hpp>

class AgradFwdLdexp : public testing::Test {
  void SetUp() {}
};

TEST_F(AgradFwdLdexp, Fvar) {
  using stan::math::exp2;
  using stan::math::fvar;
  using stan::math::ldexp;

  fvar<double> x(0.5, 1.0);

  fvar<double> a = ldexp(x, 5);
  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_);
  EXPECT_FLOAT_EQ(exp2(5), a.d_);

  fvar<double> b = 2 * ldexp(x, 5) + 4;
  EXPECT_FLOAT_EQ(2 * ldexp(0.5, 5) + 4, b.val_);
  EXPECT_FLOAT_EQ(2 * exp2(5), b.d_);

  fvar<double> c = -ldexp(x, 5) + 5;
  EXPECT_FLOAT_EQ(-ldexp(0.5, 5) + 5, c.val_);
  EXPECT_FLOAT_EQ(-exp2(5), c.d_);

  fvar<double> d = -3 * ldexp(-x, 5) + 5 * x;
  EXPECT_FLOAT_EQ(-3 * ldexp(-0.5, 5) + 5 * 0.5, d.val_);
  EXPECT_FLOAT_EQ(3 * exp2(5) + 5, d.d_);

  fvar<double> y(-0.5, 1.0);
  fvar<double> e = ldexp(y, 5);
  EXPECT_FLOAT_EQ(ldexp(-0.5, 5), e.val_);
  EXPECT_FLOAT_EQ(exp2(5), e.d_);

  fvar<double> z(0.0, 1.0);
  fvar<double> f = ldexp(z, 5);
  EXPECT_FLOAT_EQ(ldexp(0.0, 5), f.val_);
  EXPECT_FLOAT_EQ(exp2(5), f.d_);
}

TEST_F(AgradFwdLdexp, FvarFvarDouble) {
  using stan::math::exp2;
  using stan::math::fvar;
  using stan::math::ldexp;

  fvar<fvar<double> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<double> > a = ldexp(x, 5);

  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_.val_);
  EXPECT_FLOAT_EQ(exp2(5), a.val_.d_);
  EXPECT_FLOAT_EQ(0, a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);

  fvar<fvar<double> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  a = ldexp(y, 5);
  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_.val_);
  EXPECT_FLOAT_EQ(0, a.val_.d_);
  EXPECT_FLOAT_EQ(exp2(5), a.d_.val_);
  EXPECT_FLOAT_EQ(0, a.d_.d_);
}
