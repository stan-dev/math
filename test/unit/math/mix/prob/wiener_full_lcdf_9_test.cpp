#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_t0_v_sw, wiener_full_lcdf) {
  auto f_t0_v_sw = [](const auto& y, const auto& a, const auto& w,
                      const auto& sv, const auto& st0) {
    return
        [&y, &a, &w, &sv, &st0](const auto& t0, const auto& v, const auto& sw) {
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

  stan::test::expect_ad(f_t0_v_sw(y, a, w, sv, st0), t0, v, sw);
}

TEST(mathMixScalFun_t0_v_st0, wiener_full_lcdf) {
  auto f_t0_v_st0 = [](const auto& y, const auto& a, const auto& w,
                       const auto& sv, const auto& sw) {
    return
        [&y, &a, &w, &sv, &sw](const auto& t0, const auto& v, const auto& st0) {
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

  stan::test::expect_ad(f_t0_v_st0(y, a, w, sv, sw), t0, v, st0);
}

TEST(mathMixScalFun_t0_sv_sw, wiener_full_lcdf) {
  auto f_t0_sv_sw = [](const auto& y, const auto& a, const auto& w,
                       const auto& v, const auto& st0) {
    return
        [&y, &a, &w, &v, &st0](const auto& t0, const auto& sv, const auto& sw) {
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

  stan::test::expect_ad(f_t0_sv_sw(y, a, w, v, st0), t0, sv, sw);
}

TEST(mathMixScalFun_t0_sv_st0, wiener_full_lcdf) {
  auto f_t0_sv_st0 = [](const auto& y, const auto& a, const auto& w,
                        const auto& v, const auto& sw) {
    return
        [&y, &a, &w, &v, &sw](const auto& t0, const auto& sv, const auto& st0) {
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

  stan::test::expect_ad(f_t0_sv_st0(y, a, w, v, sw), t0, sv, st0);
}

TEST(mathMixScalFun_t0_sw_st0, wiener_full_lcdf) {
  auto f_t0_sw_st0 = [](const auto& y, const auto& a, const auto& w,
                        const auto& v, const auto& sv) {
    return
        [&y, &a, &w, &v, &sv](const auto& t0, const auto& sw, const auto& st0) {
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

  stan::test::expect_ad(f_t0_sw_st0(y, a, w, v, sv), t0, sw, st0);
}
