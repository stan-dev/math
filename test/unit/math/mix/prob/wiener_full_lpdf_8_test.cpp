#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_t0_w_v, wiener_full_lpdf) {
  auto f_t0_w_v = [](const auto& y, const auto& a, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&y, &a, &sv, &sw, &st0](const auto& t0, const auto& w, const auto& v) {
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

  stan::test::expect_ad(f_t0_w_v(y, a, sv, sw, st0), t0, w, v);
}

TEST(mathMixScalFun_t0_w_sv, wiener_full_lpdf) {
  auto f_t0_w_sv = [](const auto& y, const auto& a, const auto& v,
                      const auto& sw, const auto& st0) {
    return
        [&y, &a, &v, &sw, &st0](const auto& t0, const auto& w, const auto& sv) {
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

  stan::test::expect_ad(f_t0_w_sv(y, a, v, sw, st0), t0, w, sv);
}

TEST(mathMixScalFun_t0_w_sw, wiener_full_lpdf) {
  auto f_t0_w_sw = [](const auto& y, const auto& a, const auto& v,
                      const auto& sv, const auto& st0) {
    return
        [&y, &a, &v, &sv, &st0](const auto& t0, const auto& w, const auto& sw) {
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

  stan::test::expect_ad(f_t0_w_sw(y, a, v, sv, st0), t0, w, sw);
}

TEST(mathMixScalFun_t0_w_st0, wiener_full_lpdf) {
  auto f_t0_w_st0 = [](const auto& y, const auto& a, const auto& v,
                       const auto& sv, const auto& sw) {
    return
        [&y, &a, &v, &sv, &sw](const auto& t0, const auto& w, const auto& st0) {
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

  stan::test::expect_ad(f_t0_w_st0(y, a, v, sv, sw), t0, w, st0);
}

TEST(mathMixScalFun_t0_v_sv, wiener_full_lpdf) {
  auto f_t0_v_sv = [](const auto& y, const auto& a, const auto& w,
                      const auto& sw, const auto& st0) {
    return
        [&y, &a, &w, &sw, &st0](const auto& t0, const auto& v, const auto& sv) {
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

  stan::test::expect_ad(f_t0_v_sv(y, a, w, sw, st0), t0, v, sv);
}
