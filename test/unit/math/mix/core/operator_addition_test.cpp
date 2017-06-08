#include <stan/math.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>
#include <gtest/gtest.h>
#include <iostream>

#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

struct op_add_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply(const T1& x1, const T2& x2) {
    return x1 + x2;
  }
};

TEST(foo, framework) {
  using stan::math::test::test_ad;

  double nan = stan::math::not_a_number();
  double pos_inf = stan::math::positive_infinity();
  double neg_inf = stan::math::negative_infinity();

  test_ad<op_add_f>(0.5, 1.3, 1.8);
  test_ad<op_add_f>(1, 2.0, 3);
  test_ad<op_add_f>(-1.5, 1.5, 0);

  test_ad<op_add_f>(pos_inf, 0, pos_inf);
  test_ad<op_add_f>(0, pos_inf, pos_inf);
  test_ad<op_add_f>(pos_inf, pos_inf, pos_inf);

  test_ad<op_add_f>(neg_inf, 0, neg_inf);
  test_ad<op_add_f>(0, neg_inf, neg_inf);
  test_ad<op_add_f>(neg_inf, neg_inf, neg_inf);

  test_ad<op_add_f>(nan, 1.5, nan);
  test_ad<op_add_f>(0.2, nan, nan);
  test_ad<op_add_f>(nan, nan, nan);

  test_ad<op_add_f>(pos_inf, neg_inf, nan);
}




// TEST(AgradMixOperatorAddition,FvarVar_FvarVar_1stDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<var> x(0.5,1.3);
//   fvar<var> z(0.5,1.3);
//   fvar<var> a = x + z;

//   EXPECT_FLOAT_EQ(1.0, a.val_.val());
//   EXPECT_FLOAT_EQ(2.6, a.d_.val());

//   AVEC y = createAVEC(x.val_,z.val_);
//   VEC g;
//   a.val_.grad(y,g);
//   EXPECT_FLOAT_EQ(1, g[0]);
//   EXPECT_FLOAT_EQ(1, g[1]);
// }

// TEST(AgradMixOperatorAddition,FvarVar_Double_1stDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<var> x(0.5,1.3);
//   double z(0.5);
//   fvar<var> a = x + z;

//   EXPECT_FLOAT_EQ(1.0, a.val_.val());
//   EXPECT_FLOAT_EQ(1.3, a.d_.val());

//   AVEC y = createAVEC(x.val_);
//   VEC g;
//   a.val_.grad(y,g);
//   EXPECT_FLOAT_EQ(1, g[0]);
// }

// TEST(AgradMixOperatorAddition,Double_FvarVar_1stDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   double x(0.5);
//   fvar<var> z(0.5,1.3);
//   fvar<var> a = x + z;

//   EXPECT_FLOAT_EQ(1.0, a.val_.val());
//   EXPECT_FLOAT_EQ(1.3, a.d_.val());

//   AVEC y = createAVEC(z.val_);
//   VEC g;
//   a.val_.grad(y,g);
//   EXPECT_FLOAT_EQ(1, g[0]);
// }
// TEST(AgradMixOperatorAddition,FvarVar_FvarVar_2ndDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<var> x(0.5,1.3);
//   fvar<var> z(0.5,1.3);
//   fvar<var> a = x + z;

//   AVEC y = createAVEC(x.val_,z.val_);
//   VEC g;
//   a.d_.grad(y,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
//   EXPECT_FLOAT_EQ(0, g[1]);
// }
// TEST(AgradMixOperatorAddition,FvarVar_Double_2ndDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<var> x(0.5,1.3);
//   double z(0.5);
//   fvar<var> a = x + z;

//   AVEC y = createAVEC(x.val_);
//   VEC g;
//   a.d_.grad(y,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
// }
// TEST(AgradMixOperatorAddition,Double_FvarVar_2ndDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   double x(0.5);
//   fvar<var> z(0.5,1.3);
//   fvar<var> a = x + z;

//   AVEC y = createAVEC(z.val_);
//   VEC g;
//   a.d_.grad(y,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
// }
// TEST(AgradMixOperatorAddition,FvarFvarVar_FvarFvarVar_1stDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<fvar<var> > x;
//   x.val_.val_ = 0.5;
//   x.val_.d_ = 1.0;

//   fvar<fvar<var> > y;
//   y.val_.val_ = 0.5;
//   y.d_.val_ = 1.0;

//   fvar<fvar<var> > z = x + y;
//   EXPECT_FLOAT_EQ(1, z.val_.val_.val());
//   EXPECT_FLOAT_EQ(1, z.val_.d_.val());
//   EXPECT_FLOAT_EQ(1, z.d_.val_.val());
//   EXPECT_FLOAT_EQ(0, z.d_.d_.val());

//   AVEC p = createAVEC(x.val_.val_,y.val_.val_);
//   VEC g;
//   z.val_.val_.grad(p,g);
//   EXPECT_FLOAT_EQ(1, g[0]);
//   EXPECT_FLOAT_EQ(1, g[1]);
// }
// TEST(AgradMixOperatorAddition,FvarFvarVar_Double_1stDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<fvar<var> > x;
//   x.val_.val_ = 0.5;
//   x.val_.d_ = 1.0;
//   double y(0.5);

//   fvar<fvar<var> > z = x + y;
//   EXPECT_FLOAT_EQ(1, z.val_.val_.val());
//   EXPECT_FLOAT_EQ(1, z.val_.d_.val());
//   EXPECT_FLOAT_EQ(0, z.d_.val_.val());
//   EXPECT_FLOAT_EQ(0, z.d_.d_.val());

//   AVEC p = createAVEC(x.val_.val_);
//   VEC g;
//   z.val_.val_.grad(p,g);
//   EXPECT_FLOAT_EQ(1, g[0]);
// }

// TEST(AgradMixOperatorAddition,Double_FvarFvarVar_1stDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   double x(0.5);
//   fvar<fvar<var> > y;
//   y.val_.val_ = 0.5;
//   y.d_.val_ = 1.0;

//   fvar<fvar<var> > z = x + y;
//   EXPECT_FLOAT_EQ(1, z.val_.val_.val());
//   EXPECT_FLOAT_EQ(0, z.val_.d_.val());
//   EXPECT_FLOAT_EQ(1, z.d_.val_.val());
//   EXPECT_FLOAT_EQ(0, z.d_.d_.val());

//   AVEC p = createAVEC(y.val_.val_);
//   VEC g;
//   z.val_.val_.grad(p,g);
//   EXPECT_FLOAT_EQ(1, g[0]);
// }
// TEST(AgradMixOperatorAddition,FvarFvarVar_FvarFvarVar_2ndDeriv_x) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<fvar<var> > x;
//   x.val_.val_ = 0.5;
//   x.val_.d_ = 1.0;
//   fvar<fvar<var> > y;
//   y.val_.val_ = 0.5;
//   y.d_.val_ = 1.0;

//   fvar<fvar<var> > z = x + y;

//   AVEC p = createAVEC(x.val_.val_,y.val_.val_);
//   VEC g;
//   z.val_.d_.grad(p,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
//   EXPECT_FLOAT_EQ(0, g[1]);
// }
// TEST(AgradMixOperatorAddition,FvarFvarVar_FvarFvarVar_2ndDeriv_y) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<fvar<var> > x;
//   x.val_.val_ = 0.5;
//   x.val_.d_ = 1.0;
//   fvar<fvar<var> > y;
//   y.val_.val_ = 0.5;
//   y.d_.val_ = 1.0;

//   fvar<fvar<var> > z = x + y;

//   AVEC p = createAVEC(x.val_.val_,y.val_.val_);
//   VEC g;
//   z.d_.val_.grad(p,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
//   EXPECT_FLOAT_EQ(0, g[1]);
// }
// TEST(AgradMixOperatorAddition,FvarFvarVar_Double_2ndDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<fvar<var> > x;
//   x.val_.val_ = 0.5;
//   x.val_.d_ = 1.0;
//   double y(0.5);

//   fvar<fvar<var> > z = x + y;

//   AVEC p = createAVEC(x.val_.val_);
//   VEC g;
//   z.val_.d_.grad(p,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
// }

// TEST(AgradMixOperatorAddition,Double_FvarFvarVar_2ndDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   double x(0.5);
//   fvar<fvar<var> > y;
//   y.val_.val_ = 0.5;
//   y.d_.val_ = 1.0;

//   fvar<fvar<var> > z = x + y;

//   AVEC p = createAVEC(y.val_.val_);
//   VEC g;
//   z.d_.val_.grad(p,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
// }
// TEST(AgradMixOperatorAddition,FvarFvarVar_FvarFvarVar_3rdDeriv_y) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<fvar<var> > x;
//   x.val_.val_ = 0.5;
//   x.val_.d_ = 1.0;
//   x.d_.val_ = 1.0;
//   fvar<fvar<var> > y;
//   y.val_.val_ = 0.5;
//   y.d_.val_ = 1.0;
//   y.val_.d_ = 1.0;

//   fvar<fvar<var> > z = x + y;

//   AVEC p = createAVEC(x.val_.val_,y.val_.val_);
//   VEC g;
//   z.d_.d_.grad(p,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
//   EXPECT_FLOAT_EQ(0, g[1]);
// }
// TEST(AgradMixOperatorAddition,FvarFvarVar_Double_3rdDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   fvar<fvar<var> > x;
//   x.val_.val_ = 0.5;
//   x.val_.d_ = 1.0;
//   x.d_.val_ = 1.0;
//   double y(0.5);

//   fvar<fvar<var> > z = x + y;

//   AVEC p = createAVEC(x.val_.val_);
//   VEC g;
//   z.d_.d_.grad(p,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
// }

// TEST(AgradMixOperatorAddition,Double_FvarFvarVar_3rdDeriv) {
//   using stan::math::fvar;
//   using stan::math::var;

//   double x(0.5);
//   fvar<fvar<var> > y;
//   y.val_.val_ = 0.5;
//   y.d_.val_ = 1.0;
//   y.val_.d_ = 1.0;

//   fvar<fvar<var> > z = x + y;

//   AVEC p = createAVEC(y.val_.val_);
//   VEC g;
//   z.d_.d_.grad(p,g);
//   EXPECT_FLOAT_EQ(0, g[0]);
// }

// struct add_fun {
//   template <typename T0, typename T1>
//   inline
//   typename stan::return_type<T0,T1>::type
//   operator()(const T0& arg1,
//              const T1& arg2) const {
//     return arg1+arg2;
//   }
// };

// TEST(AgradMixOperatorAddition, add_nan) {
//   add_fun add_;
//   test_nan_mix(add_,3.0,5.0,false);
// }
