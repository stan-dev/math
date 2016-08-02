#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

   

TEST(AgradFwdInvSquare,FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_square;

  fvar<var> x(0.5,1.0);
  fvar<var> a = inv_square(x);

  EXPECT_FLOAT_EQ(inv_square(0.5), a.val_.val());
  EXPECT_FLOAT_EQ(-2.0 * inv_square(0.5) / (0.5), a.d_.val());

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.val_.grad(y,g);
  EXPECT_FLOAT_EQ(-2.0 / (0.5 * 0.5 * 0.5), g[0]);
}
TEST(AgradFwdInvSquare,FvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_square;

  fvar<var> x(0.5,1.0);
  fvar<var> a = inv_square(x);

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.d_.grad(y,g);
  EXPECT_FLOAT_EQ(-2.0 * -3.0 / (0.5 * 0.5 * 0.5 * 0.5), g[0]);
}


TEST(AgradFwdInvSquare,FvarFvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_square;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = inv_square(x);

  EXPECT_FLOAT_EQ(inv_square(0.5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(-2.0 * inv_square(0.5) / (0.5), a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(-2.0 * inv_square(0.5) / (0.5), g[0]);
}
TEST(AgradFwdInvSquare,FvarFvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inv_square;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = inv_square(x);

  EXPECT_FLOAT_EQ(inv_square(0.5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(-2.0 * inv_square(0.5) / (0.5), a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(-2.0 * -3.0 / (0.5 * 0.5 * 0.5 * 0.5), g[0]);
}
TEST(AgradFwdInvSquare,FvarFvarVar_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  x.d_.val_ = 1.0;

  fvar<fvar<var> > a = inv_square(x);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.d_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(-768, g[0]);
}

struct inv_square_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return inv_square(arg1);
  }
};

TEST(AgradFwdInvSquare,inv_square_NaN) {
  inv_square_fun inv_square_;
  test_nan_mix(inv_square_,false);
}
