#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>

TEST(AgradFwdOwensT, FvarVar_FvarVar_FvarVar_1stDeriv) {
  using boost::math::owens_t;
  using stan::math::fvar;
  using stan::math::binormal_integral_owens;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(1.0, 1.0);
  fvar<var> rho(0.7, 1.0);
  fvar<var> f = binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.8370300138577811, f.val_.val());

}
TEST(AgradFwdOwensT, Double_FvarVar_FvarVar_1stDeriv) {
  using boost::math::owens_t;
  using stan::math::fvar;
  using stan::math::binormal_integral_owens;
  using stan::math::var;

  double a(2.0);
  fvar<var> b(1.0, 1.0);
  fvar<var> rho(0.7, 1.0);
  fvar<var> f = binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.8370300138577811, f.val_.val());
}
TEST(AgradFwdOwensT, FvarVar_Double_FvarVar_1stDeriv) {
  using boost::math::owens_t;
  using stan::math::fvar;
  using stan::math::binormal_integral_owens;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  double b = 1.0;
  fvar<var> rho(0.7, 1.0);
  fvar<var> f = binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.8370300138577811, f.val_.val());

}
TEST(AgradFwdOwensT, FvarVar_FvarVar_Double_1stDeriv) {
  using boost::math::owens_t;
  using stan::math::fvar;
  using stan::math::binormal_integral_owens;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  fvar<var> b(1.0, 1.0);
  double rho(0.7);
  fvar<var> f = binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.8370300138577811, f.val_.val());

}
TEST(AgradFwdOwensT, Double_FvarVar_Double_1stDeriv) {
  using boost::math::owens_t;
  using stan::math::fvar;
  using stan::math::binormal_integral_owens;
  using stan::math::var;

  double a(2.0);
  fvar<var> b(1.0, 1.0);
  double rho(0.7);
  fvar<var> f = binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.8370300138577811, f.val_.val());

}
TEST(AgradFwdOwensT, FvarVar_Double_Double_1stDeriv) {
  using boost::math::owens_t;
  using stan::math::fvar;
  using stan::math::binormal_integral_owens;
  using stan::math::var;

  fvar<var> a(2.0, 1.0);
  double b(1.0);
  double rho(0.7);
  fvar<var> f = binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.8370300138577811, f.val_.val());

}
TEST(AgradFwdOwensT, Double_Double_FvarVar_1stDeriv) {
  using boost::math::owens_t;
  using stan::math::fvar;
  using stan::math::binormal_integral_owens;
  using stan::math::var;

  double a(2.0);
  double b(1.0);
  fvar<var> rho(0.7, 1.0);
  fvar<var> f = binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.8370300138577811, f.val_.val());

}
