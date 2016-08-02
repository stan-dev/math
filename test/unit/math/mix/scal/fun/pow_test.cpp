#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

TEST(AgradFwdPow, FvarVar_FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::log;
  using std::pow;

  fvar<var> x(0.5,1.0);
  fvar<var> z(1.2,2.0);
  fvar<var> a = pow(x,z);

  EXPECT_FLOAT_EQ(pow(0.5,1.2), a.val_.val());
  EXPECT_FLOAT_EQ((2.0 * log(0.5) + 1.2 * 1.0 / 0.5) * pow(0.5, 1.2), 
                  a.d_.val());

  AVEC y = createAVEC(x.val_,z.val_);
  VEC g;
  a.val_.grad(y,g);
  EXPECT_FLOAT_EQ(1.2 / 0.5 * pow(0.5, 1.2), g[0]);
  EXPECT_FLOAT_EQ(log(0.5) * pow(0.5, 1.2), g[1]);
}
TEST(AgradFwdPow, FvarVar_Double_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::log;
  using std::pow;

  fvar<var> x(0.5,1.0);
  double z(1.2);
  fvar<var> a = pow(x,z);

  EXPECT_FLOAT_EQ(pow(0.5,1.2), a.val_.val());
  EXPECT_FLOAT_EQ((1.2 * 1.0 / 0.5) * pow(0.5, 1.2), 
                  a.d_.val());

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.val_.grad(y,g);
  EXPECT_FLOAT_EQ(1.2 / 0.5 * pow(0.5, 1.2), g[0]);
}
TEST(AgradFwdPow, Double_FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::log;
  using std::pow;

  double x(0.5);
  fvar<var> z(1.2,2.0);
  fvar<var> a = pow(x,z);

  EXPECT_FLOAT_EQ(pow(0.5,1.2), a.val_.val());
  EXPECT_FLOAT_EQ((2.0 * log(0.5)) * pow(0.5, 1.2), 
                  a.d_.val());

  AVEC y = createAVEC(z.val_);
  VEC g;
  a.val_.grad(y,g);
  EXPECT_FLOAT_EQ(log(0.5) * pow(0.5, 1.2), g[0]);
}
TEST(AgradFwdPow, FvarVar_FvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::log;
  using std::pow;

  fvar<var> x(0.5,1.0);
  fvar<var> z(1.2,1.0);
  fvar<var> a = pow(x,z);

  AVEC y = createAVEC(x.val_,z.val_);
  VEC g;
  a.d_.grad(y,g);
  EXPECT_FLOAT_EQ(0.56431121, g[0]);
  EXPECT_FLOAT_EQ(0.35557628, g[1]);
}
TEST(AgradFwdPow, FvarVar_Double_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::log;
  using std::pow;

  fvar<var> x(0.5,1.0);
  double z(1.2);
  fvar<var> a = pow(x,z);

  AVEC y = createAVEC(x.val_);
  VEC g;
  a.d_.grad(y,g);
  EXPECT_FLOAT_EQ((1.2 - 1.0) * 1.2 * pow(0.5,1.2 - 2.0), g[0]);
}
TEST(AgradFwdPow, Double_FvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::log;
  using std::pow;

  double x(0.5);
  fvar<var> z(1.2,1.0);
  fvar<var> a = pow(x,z);

  AVEC y = createAVEC(z.val_);
  VEC g;
  a.d_.grad(y,g);
  EXPECT_FLOAT_EQ(pow(0.5,1.2) * log(0.5) * log(0.5), g[0]);
}

TEST(AgradFwdPow, FvarFvarVar_FvarFvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = pow(x,y);

  EXPECT_FLOAT_EQ(pow(0.5,0.5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(0.5 * pow(0.5,-0.5), a.val_.d_.val());
  EXPECT_FLOAT_EQ(log(0.5) * pow(0.5,0.5), a.d_.val_.val());
  EXPECT_FLOAT_EQ(pow(0.5, -0.5) * (0.5 * log(0.5) + 1.0), a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_,y.val_.val_);
  VEC g;
  a.val_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(0.5 * pow(0.5,-0.5), g[0]);
  EXPECT_FLOAT_EQ(log(0.5) * pow(0.5,0.5), g[1]);
}
TEST(AgradFwdPow, FvarFvarVar_Double_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  double y(0.5);

  fvar<fvar<var> > a = pow(x,y);

  EXPECT_FLOAT_EQ(pow(0.5,0.5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(0.5 * pow(0.5,-0.5), a.val_.d_.val());
  EXPECT_FLOAT_EQ(0, a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(0.5 * pow(0.5,-0.5), g[0]);
}
TEST(AgradFwdPow, Double_FvarFvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  double x(0.5);
  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = pow(x,y);

  EXPECT_FLOAT_EQ(pow(0.5,0.5), a.val_.val_.val());
  EXPECT_FLOAT_EQ(0, a.val_.d_.val());
  EXPECT_FLOAT_EQ(log(0.5) * pow(0.5,0.5), a.d_.val_.val());
  EXPECT_FLOAT_EQ(0, a.d_.d_.val());

  AVEC p = createAVEC(y.val_.val_);
  VEC g;
  a.val_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(log(0.5) * pow(0.5,0.5), g[0]);
}
TEST(AgradFwdPow, FvarFvarVar_FvarFvarVar_2ndDeriv_x) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = pow(x,y);

  AVEC p = createAVEC(x.val_.val_,y.val_.val_);
  VEC g;
  a.val_.d_.grad(p,g);
  EXPECT_FLOAT_EQ((0.5 - 1.0) * 0.5 * pow(0.5,0.5 - 2.0), g[0]);
  EXPECT_FLOAT_EQ(pow(0.5, -0.5) * (0.5 * log(0.5) + 1.0), g[1]);
}
TEST(AgradFwdPow, FvarFvarVar_FvarFvarVar_2ndDeriv_y) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = pow(x,y);

  AVEC p = createAVEC(x.val_.val_,y.val_.val_);
  VEC g;
  a.d_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(pow(0.5, -0.5) * (0.5 * log(0.5) + 1.0), g[0]);
  EXPECT_FLOAT_EQ(pow(0.5,0.5) * log(0.5) * log(0.5), g[1]);
}
TEST(AgradFwdPow, FvarFvarVar_Double_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  double y(1.2);

  fvar<fvar<var> > a = pow(x,y);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.val_.d_.grad(p,g);
  EXPECT_FLOAT_EQ((1.2 - 1.0) * 1.2 * pow(0.5,1.2 - 2.0), g[0]);
}
TEST(AgradFwdPow, Double_FvarFvarVar_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  double x(0.5);
  fvar<fvar<var> > y;
  y.val_.val_ = 1.2;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = pow(x,y);

  AVEC p = createAVEC(y.val_.val_);
  VEC g;
  a.d_.val_.grad(p,g);
  EXPECT_FLOAT_EQ(pow(0.5,1.2) * log(0.5) * log(0.5), g[0]);
}
TEST(AgradFwdPow, FvarFvarVar_FvarFvarVar_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  fvar<fvar<var> > y;
  y.val_.val_ = 0.5;
  y.d_.val_ = 1.0;

  fvar<fvar<var> > a = pow(x,y);

  AVEC p = createAVEC(x.val_.val_,y.val_.val_);
  VEC g;
  a.d_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(0.49012908, g[0]);
  EXPECT_FLOAT_EQ(-1.6207848, g[1]);
}
TEST(AgradFwdPow, FvarFvarVar_Double_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  fvar<fvar<var> > x;
  x.val_.val_ = 0.5;
  x.val_.d_ = 1.0;
  x.d_.val_ = 1.0;
  double y(1.2);

  fvar<fvar<var> > a = pow(x,y);

  AVEC p = createAVEC(x.val_.val_);
  VEC g;
  a.d_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(-0.668583, g[0]);
}
TEST(AgradFwdPow, Double_FvarFvarVar_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using std::pow;
  using std::log;

  double x(0.5);
  fvar<fvar<var> > y;
  y.val_.val_ = 1.2;
  y.d_.val_ = 1.0;
  y.val_.d_ = 1.0;

  fvar<fvar<var> > a = pow(x,y);

  AVEC p = createAVEC(y.val_.val_);
  VEC g;
  a.d_.d_.grad(p,g);
  EXPECT_FLOAT_EQ(-0.14495739, g[0]);
}

struct pow_fun {
  template <typename T0, typename T1>
  inline 
  typename boost::math::tools::promote_args<T0,T1>::type
  operator()(const T0 arg1,
             const T1 arg2) const {
    return pow(arg1,arg2);
  }
};

TEST(AgradFwdPow, nan) {
  pow_fun pow_;
  test_nan_mix(pow_,3.0,5.0,false);
}
