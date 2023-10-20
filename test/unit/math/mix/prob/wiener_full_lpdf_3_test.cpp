#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_y_w_v, wiener_full_lpdf) {
  auto f_y_w_v = [](const auto& a, const auto& t0, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&a, &t0, &sv, &sw, &st0](const auto& y, const auto& w, const auto& v) {
          return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_y_w_v(a, t0, sv, sw, st0), y, w, v);
}

TEST(mathMixScalFun_y_w_sv, wiener_full_lpdf) {
  auto f_y_w_sv = [](const auto& a, const auto& t0, const auto& v,
                     const auto& sw, const auto& st0) {
    return
        [&a, &t0, &v, &sw, &st0](const auto& y, const auto& w, const auto& sv) {
          return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_y_w_sv(a, t0, v, sw, st0), y, w, sv);
}

TEST(mathMixScalFun_y_w_sw, wiener_full_lpdf) {
  auto f_y_w_sw = [](const auto& a, const auto& t0, const auto& v,
                     const auto& sv, const auto& st0) {
    return
        [&a, &t0, &v, &sv, &st0](const auto& y, const auto& w, const auto& sw) {
          return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_y_w_sw(a, t0, v, sv, st0), y, w, sw);
}

TEST(mathMixScalFun_y_w_st0, wiener_full_lpdf) {
  auto f_y_w_st0 = [](const auto& a, const auto& t0, const auto& v,
                      const auto& sv, const auto& sw) {
    return
        [&a, &t0, &v, &sv, &sw](const auto& y, const auto& w, const auto& st0) {
          return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_y_w_st0(a, t0, v, sv, sw), y, w, st0);
}

TEST(mathMixScalFun_y_v_sv, wiener_full_lpdf) {
  auto f_y_v_sv = [](const auto& a, const auto& t0, const auto& w,
                     const auto& sw, const auto& st0) {
    return
        [&a, &t0, &w, &sw, &st0](const auto& y, const auto& v, const auto& sv) {
          return stan::math::wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_y_v_sv(a, t0, w, sw, st0), y, v, sv);
}
