#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/erf.hpp>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>



TEST(AgradFwdErfc,FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  using boost::math::erfc;

  fvar<var> x(0.5,1.3);
  fvar<var> a = erfc(x);

  EXPECT_FLOAT_EQ(erfc(0.5), a.val_.val());
  EXPECT_FLOAT_EQ(-1.3 * 2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), a.d_.val());

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.val_.grad(y,g);
  EXPECT_FLOAT_EQ(-2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), g[0]);
}
TEST(AgradFwdErfc,FvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  using boost::math::erfc;

  fvar<var> x(0.5,1.3);
  fvar<var> a = erfc(x);

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.d_.grad(y,g);
  EXPECT_FLOAT_EQ(1.3 * 2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), g[0]);
}



TEST(AgradFwdErfc,FvarFvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  using boost::math::erfc;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = erfc(x);

  EXPECT_FLOAT_EQ(erfc(0.5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(-2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(-2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), g[0]);

  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > b = erfc(y);
  EXPECT_FLOAT_EQ(erfc(0.5), b.val_.val_.val());
  EXPECT_FLOAT_EQ(0, b.val_.d_.val());
  EXPECT_FLOAT_EQ(-2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), b.d_.val_.val());
  EXPECT_FLOAT_EQ(0, b.d_.d_.val());


  AVEC q = createAVEC(y.val_.val_);
  VEC r;
  b.val_.val_.grad(q,r);
  EXPECT_FLOAT_EQ(-2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), r[0]);
}

TEST(AgradFwdErfc,FvarFvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  using boost::math::erfc;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;

  fvar<fvar<var> > a = erfc(x);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), g[0]);

  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > b = erfc(y);

  AVEC q = createAVEC(y.val_.val_);
  VEC r;
  b.d_.val_.grad(q,r);
  EXPECT_FLOAT_EQ(2 * exp(-0.5 * 0.5) / 
                  sqrt(boost::math::constants::pi<double>()), r[0]);
}
TEST(AgradFwdErfc,FvarFvarVar_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::exp;
  using std::sqrt;
  using boost::math::erfc;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  x.d_.val_ = 1.0;

  fvar<fvar<var> > a = erfc(x);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.d_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(0.878782578935444794093723954824, g[0]);
}

struct erfc_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return erfc(arg1);
  }
};

TEST(AgradFwdErfc,erfc_NaN) {
  erfc_fun erfc_;
  test_nan_mix(erfc_,false);
}
