#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_a_v_sw, wiener_full_lcdf) {
  auto f_a_v_sw = [](const auto& y, const auto& t0, const auto& w,
                     const auto& sv, const auto& st0) {
    return
        [&y, &t0, &w, &sv, &st0](const auto& a, const auto& v, const auto& sw) {
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

  stan::test::expect_ad(f_a_v_sw(y, t0, w, sv, st0), a, v, sw);
}

TEST(mathMixScalFun_a_v_st0, wiener_full_lcdf) {
  auto f_a_v_st0 = [](const auto& y, const auto& t0, const auto& w,
                      const auto& sv, const auto& sw) {
    return
        [&y, &t0, &w, &sv, &sw](const auto& a, const auto& v, const auto& st0) {
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

  stan::test::expect_ad(f_a_v_st0(y, t0, w, sv, sw), a, v, st0);
}

TEST(mathMixScalFun_a_sv_sw, wiener_full_lcdf) {
  auto f_a_sv_sw = [](const auto& y, const auto& t0, const auto& w,
                      const auto& v, const auto& st0) {
    return
        [&y, &t0, &w, &v, &st0](const auto& a, const auto& sv, const auto& sw) {
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

  stan::test::expect_ad(f_a_sv_sw(y, t0, w, v, st0), a, sv, sw);
}

TEST(mathMixScalFun_a_sv_st0, wiener_full_lcdf) {
  auto f_a_sv_st0 = [](const auto& y, const auto& t0, const auto& w,
                       const auto& v, const auto& sw) {
    return
        [&y, &t0, &w, &v, &sw](const auto& a, const auto& sv, const auto& st0) {
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

  stan::test::expect_ad(f_a_sv_st0(y, t0, w, v, sw), a, sv, st0);
}

TEST(mathMixScalFun_a_sw_st0, wiener_full_lcdf) {
  auto f_a_sw_st0 = [](const auto& y, const auto& t0, const auto& w,
                       const auto& v, const auto& sv) {
    return
        [&y, &t0, &w, &v, &sv](const auto& a, const auto& sw, const auto& st0) {
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

  stan::test::expect_ad(f_a_sw_st0(y, t0, w, v, sv), a, sw, st0);
}
