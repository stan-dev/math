#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>


TEST(AgradFwdInvLogit,FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_logit;

  fvar<var> x(0.5,1.3);
  fvar<var> a = inv_logit(x);

  EXPECT_FLOAT_EQ(inv_logit(0.5), a.val_.val());
  EXPECT_FLOAT_EQ(1.3 * inv_logit(0.5) * (1 - inv_logit(0.5)), a.d_.val());

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.val_.grad(y,g);
  EXPECT_FLOAT_EQ(inv_logit(0.5) * (1 - inv_logit(0.5)), g[0]);
}
TEST(AgradFwdInvLogit,FvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_logit;

  fvar<var> x(0.5,1.3);
  fvar<var> a = inv_logit(x);

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.d_.grad(y,g);
  EXPECT_FLOAT_EQ(1.3 * (inv_logit(0.5) * (1 - inv_logit(0.5)) - inv_logit(0.5) 
                         * 2.0 * inv_logit(0.5) * (1 - inv_logit(0.5))), g[0]);
}

TEST(AgradFwdInvLogit,FvarFvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_logit;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = inv_logit(x);

  EXPECT_FLOAT_EQ(inv_logit(0.5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(inv_logit(0.5) * (1 - inv_logit(0.5)), a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(inv_logit(0.5) * (1 - inv_logit(0.5)), g[0]);

  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > b = inv_logit(y);
  EXPECT_FLOAT_EQ(inv_logit(0.5), b.val_.val_.val());
  EXPECT_FLOAT_EQ(0, b.val_.d_.val());
  EXPECT_FLOAT_EQ(inv_logit(0.5) * (1 - inv_logit(0.5)), b.d_.val_.val());
  EXPECT_FLOAT_EQ(0, b.d_.d_.val());

  AVEC q = createAVEC(y.val_.val_);
  VEC r;
  b.val_.val_.grad(q,r);
  EXPECT_FLOAT_EQ(inv_logit(0.5) * (1 - inv_logit(0.5)), r[0]);
}
TEST(AgradFwdInvLogit,FvarFvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_logit;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = inv_logit(x);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.d_.grad(p,g);
  EXPECT_FLOAT_EQ((inv_logit(0.5) * (1 - inv_logit(0.5)) - inv_logit(0.5) 
                   * 2.0 * inv_logit(0.5) * (1 - inv_logit(0.5))), g[0]);

  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > b = inv_logit(y);

  AVEC q = createAVEC(y.val_.val_);
  VEC r;
  b.d_.val_.grad(q,r);
  EXPECT_FLOAT_EQ((inv_logit(0.5) * (1 - inv_logit(0.5)) - inv_logit(0.5) 
                   * 2.0 * inv_logit(0.5) * (1 - inv_logit(0.5))), r[0]);
}
TEST(AgradFwdInvLogit,FvarFvarVar_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_logit;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  x.d_.val_ = 1.0;

  fvar<fvar<var> > a = inv_logit(x);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.d_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(-0.09635675628958461417489403610, g[0]);
}

struct inv_logit_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return inv_logit(arg1);
  }
};

TEST(AgradFwdInvLogit,inv_logit_NaN) {
  inv_logit_fun inv_logit_;
  test_nan_mix(inv_logit_,false);
}
