#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixDouble5, wiener_lpdf) {
  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 1.5;
  double sv = 0.2;
  stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
}

TEST(mathMixVar5, wiener_lpdf) {
  using stan::math::var;
  var y = 1.0;
  var a = 2.0;
  var t0 = 0.2;
  var w = 0.5;
  var v = 1.5;
  var sv = 0.2;
  stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
}

TEST(mathMixFVar5, wiener_lpdf) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> y = 1.0;
  fvar<var> a = 2.0;
  fvar<var> t0 = 0.2;
  fvar<var> w = 0.5;
  fvar<var> v = 1.5;
  fvar<var> sv = 0.2;
  double error = 1e-4;
  stan::math::wiener_lpdf(y, a, t0, w, v, sv, error);
}
