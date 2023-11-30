#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixDouble, wiener_full_lcdf) {
  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 1.5;
  double sv = 0.2;
  double sw = 0.2;
  double st0 = 0.2;
  stan::math::wiener_full_lcdf(y, a, t0, w, v, sv, sw, st0, 1e-4);
}

TEST(mathMixVar, wiener_full_lcdf) {
  using stan::math::var;
  var y = 1.0;
  var a = 2.0;
  var t0 = 0.2;
  var w = 0.5;
  var v = 1.5;
  var sv = 0.2;
  var sw = 0;
  var st0 = 0.2;
  stan::math::wiener_full_lcdf(y, a, t0, w, v, sv, sw, st0);
}

TEST(mathMixFVar, wiener_full_lcdf) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> y = 1.0;
  fvar<var> a = 2.0;
  fvar<var> t0 = 0.2;
  fvar<var> w = 0.5;
  fvar<var> v = 1.5;
  fvar<var> sv = 0.2;
  fvar<var> sw = 0;
  fvar<var> st0 = 0.2;
  stan::math::wiener_full_lcdf(y, a, t0, w, v, sv, sw, st0);
}

TEST(mathMixScalFun_y_a_t0, wiener_full_lcdf) {
  auto f_y_a_t0 = [](const auto& w, const auto& v, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&w, &v, &sv, &sw, &st0](const auto& y, const auto& a, const auto& t0) {
          return stan::math::wiener_full_lcdf(y, a, t0, w, v, sv, sw, st0);
        };
  };

  double y = 0.1;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;
  double sw = 0.1;
  double st0 = 0.3;

  stan::test::expect_ad(f_y_a_t0(w, v, sv, sw, st0), y, a, t0);
}