#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

class AgradFwdLdexp : public testing::Test {
  void SetUp() { stan::math::recover_memory(); }
};

TEST_F(AgradFwdLdexp, FvarVar_1stDeriv) {
  using stan::math::exp2;
  using stan::math::fvar;
  using stan::math::ldexp;
  using stan::math::var;

  fvar<var> x(0.5, 1.3);
  fvar<var> a = ldexp(x, 5);

  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_.val());
  EXPECT_FLOAT_EQ(1.3 * exp2(5), a.d_.val());

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.val_.grad(y, g);
  EXPECT_FLOAT_EQ(exp2(5), g[0]);
}

TEST_F(AgradFwdLdexp, FvarVar_2ndDeriv) {
  using stan::math::exp2;
  using stan::math::fvar;
  using stan::math::ldexp;
  using stan::math::var;

  fvar<var> x(0.5, 1.3);
  fvar<var> a = ldexp(x, 5);

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.d_.grad(y, g);
  EXPECT_FLOAT_EQ(0.0, g[0]);
}

TEST_F(AgradFwdLdexp, FvarFvarVar_1stDeriv) {
  using stan::math::exp2;
  using stan::math::fvar;
  using stan::math::ldexp;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = ldexp(x, 5);

  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(exp2(5), a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.val_.grad(p, g);
  stan::math::recover_memory();
  EXPECT_FLOAT_EQ(exp2(5), g[0]);

  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > b = ldexp(y, 5);
  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(0, a.val_.d_.val());
  EXPECT_FLOAT_EQ(exp2(5), a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC q = createAVEC(y.val_.val_);
  VEC r;
  b.val_.val_.grad(q, r);
  EXPECT_FLOAT_EQ(exp2(5), r[0]);
}

TEST_F(AgradFwdLdexp, FvarFvarVar_2ndDeriv) {
  using stan::math::exp2;
  using stan::math::fvar;
  using stan::math::ldexp;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = ldexp(x, 5);

  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(exp2(5), a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.d_.grad(p, g);
  stan::math::recover_memory();

  EXPECT_FLOAT_EQ(0.0, g[0]);

  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > b = ldexp(y, 5);
  EXPECT_FLOAT_EQ(ldexp(0.5, 5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(0, a.val_.d_.val());
  EXPECT_FLOAT_EQ(exp2(5), a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC q = createAVEC(y.val_.val_);
  VEC r;
  b.d_.val_.grad(q, r);
  EXPECT_FLOAT_EQ(0.0, r[0]);
}
TEST_F(AgradFwdLdexp, FvarFvarVar_3rdDeriv) {
  using stan::math::exp2;
  using stan::math::fvar;
  using stan::math::ldexp;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  x.d_.val_ = 1.0;

  fvar<fvar<var> > a = ldexp(x, 5);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.d_.d_.grad(p, g);
  EXPECT_FLOAT_EQ(0.0, g[0]);
}
