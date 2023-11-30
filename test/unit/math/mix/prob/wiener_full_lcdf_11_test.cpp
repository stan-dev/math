#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_w_sw_st0, wiener_full_lcdf) {
  auto f_w_sw_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& v, const auto& sv) {
    return
        [&y, &a, &t0, &v, &sv](const auto& w, const auto& sw, const auto& st0) {
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

  stan::test::expect_ad(f_w_sw_st0(y, a, t0, v, sv), w, sw, st0);
}

TEST(mathMixScalFun_v_sv_sw, wiener_full_lcdf) {
  auto f_v_sv_sw = [](const auto& y, const auto& a, const auto& t0,
                      const auto& w, const auto& st0) {
    return
        [&y, &a, &t0, &w, &st0](const auto& v, const auto& sv, const auto& sw) {
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

  stan::test::expect_ad(f_v_sv_sw(y, a, t0, w, st0), v, sv, sw);
}

TEST(mathMixScalFun_v_sv_st0, wiener_full_lcdf) {
  auto f_v_sv_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& w, const auto& sw) {
    return
        [&y, &a, &t0, &w, &sw](const auto& v, const auto& sv, const auto& st0) {
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

  stan::test::expect_ad(f_v_sv_st0(y, a, t0, w, sw), v, sv, st0);
}

TEST(mathMixScalFun_v_sw_st0, wiener_full_lcdf) {
  auto f_v_sw_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& w, const auto& sv) {
    return
        [&y, &a, &t0, &w, &sv](const auto& v, const auto& sw, const auto& st0) {
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

  stan::test::expect_ad(f_v_sw_st0(y, a, t0, w, sv), v, sw, st0);
}

TEST(mathMixScalFun_sv_sw_st0, wiener_full_lcdf) {
  auto f_sv_sw_st0 = [](const auto& y, const auto& a, const auto& t0,
                        const auto& w, const auto& v) {
    return
        [&y, &a, &t0, &w, &v](const auto& sv, const auto& sw, const auto& st0) {
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

  stan::test::expect_ad(f_sv_sw_st0(y, a, t0, w, v), sv, sw, st0);
}
