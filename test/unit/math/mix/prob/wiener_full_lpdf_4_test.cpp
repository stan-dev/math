#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_y_v_sw, wiener_full_lpdf) {
  auto f_y_v_sw = [](const auto& a, const auto& t0, const auto& w,
                     const auto& sv, const auto& st0) {
    return
        [&a, &t0, &w, &sv, &st0](const auto& y, const auto& v, const auto& sw) {
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

  stan::test::expect_ad(f_y_v_sw(a, t0, w, sv, st0), y, v, sw);
}

TEST(mathMixScalFun_y_v_st0, wiener_full_lpdf) {
  auto f_y_v_st0 = [](const auto& a, const auto& t0, const auto& w,
                      const auto& sv, const auto& sw) {
    return
        [&a, &t0, &w, &sv, &sw](const auto& y, const auto& v, const auto& st0) {
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

  stan::test::expect_ad(f_y_v_st0(a, t0, w, sv, sw), y, v, st0);
}

TEST(mathMixScalFun_y_sv_sw, wiener_full_lpdf) {
  auto f_y_sv_sw = [](const auto& a, const auto& t0, const auto& w,
                      const auto& v, const auto& st0) {
    return
        [&a, &t0, &w, &v, &st0](const auto& y, const auto& sv, const auto& sw) {
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

  stan::test::expect_ad(f_y_sv_sw(a, t0, w, v, st0), y, sv, sw);
}

TEST(mathMixScalFun_y_sv_st0, wiener_full_lpdf) {
  auto f_y_sv_st0 = [](const auto& a, const auto& t0, const auto& w,
                       const auto& v, const auto& sw) {
    return
        [&a, &t0, &w, &v, &sw](const auto& y, const auto& sv, const auto& st0) {
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

  stan::test::expect_ad(f_y_sv_st0(a, t0, w, v, sw), y, sv, st0);
}

TEST(mathMixScalFun_y_sw_st0, wiener_full_lpdf) {
  auto f_y_sw_st0 = [](const auto& a, const auto& t0, const auto& w,
                       const auto& v, const auto& sv) {
    return
        [&a, &t0, &w, &v, &sv](const auto& y, const auto& sw, const auto& st0) {
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

  stan::test::expect_ad(f_y_sw_st0(a, t0, w, v, sv), y, sw, st0);
}
