#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixDouble, wiener4_lccdf) {
  double y = 1.0;
  double a = 2.5;
  double t0 = 0.2;
  double w = 0.5;
  double v = 1.5;
  stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
}

TEST(mathMixVar, wiener4_lccdf) {
  using stan::math::var;
  var y = 1.0;
  var a = 2.5;
  var t0 = 0.2;
  var w = 0.5;
  var v = 1.5;
  stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
}

TEST(mathMixFVar, wiener4_lccdf) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> y = 1.0;
  fvar<var> a = 2.5;
  fvar<var> t0 = 0.2;
  fvar<var> w = 0.5;
  fvar<var> v = 1.5;
  double error = 1e-4;
  stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, error);
}
