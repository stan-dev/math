#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_w_v_sv, wiener_full_lcdf) {
  auto f_w_v_sv = [](const auto& y, const auto& a, const auto& t0,
                     const auto& sw, const auto& st0) {
    return
        [&y, &a, &t0, &sw, &st0](const auto& w, const auto& v, const auto& sv) {
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

  stan::test::expect_ad(f_w_v_sv(y, a, t0, sw, st0), w, v, sv);
}

TEST(mathMixScalFun_w_v_sw, wiener_full_lcdf) {
  auto f_w_v_sw = [](const auto& y, const auto& a, const auto& t0,
                     const auto& sv, const auto& st0) {
    return
        [&y, &a, &t0, &sv, &st0](const auto& w, const auto& v, const auto& sw) {
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

  stan::test::expect_ad(f_w_v_sw(y, a, t0, sv, st0), w, v, sw);
}

TEST(mathMixScalFun_w_v_st0, wiener_full_lcdf) {
  auto f_w_v_st0 = [](const auto& y, const auto& a, const auto& t0,
                      const auto& sv, const auto& sw) {
    return
        [&y, &a, &t0, &sv, &sw](const auto& w, const auto& v, const auto& st0) {
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

  stan::test::expect_ad(f_w_v_st0(y, a, t0, sv, sw), w, v, st0);
}

TEST(mathMixScalFun_w_sv_sw, wiener_full_lcdf) {
  auto f_w_sv_sw = [](const auto& y, const auto& a, const auto& t0,
                      const auto& v, const auto& st0) {
    return
        [&y, &a, &t0, &v, &st0](const auto& w, const auto& sv, const auto& sw) {
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

  stan::test::expect_ad(f_w_sv_sw(y, a, t0, v, st0), w, sv, sw);
}

TEST(mathMixScalFun_w_sv_st0, wiener_full_lcdf) {
  auto f_w_sv_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& v, const auto& sw) {
    return
        [&y, &a, &t0, &v, &sw](const auto& w, const auto& sv, const auto& st0) {
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

  stan::test::expect_ad(f_w_sv_st0(y, a, t0, v, sw), w, sv, st0);
}
