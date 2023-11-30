#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_a_w_v, wiener_full_lcdf) {
  auto f_a_w_v = [](const auto& y, const auto& t0, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&y, &t0, &sv, &sw, &st0](const auto& a, const auto& w, const auto& v) {
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

  stan::test::expect_ad(f_a_w_v(y, t0, sv, sw, st0), a, w, v);
}

TEST(mathMixScalFun_a_w_sv, wiener_full_lcdf) {
  auto f_a_w_sv = [](const auto& y, const auto& t0, const auto& v,
                     const auto& sw, const auto& st0) {
    return
        [&y, &t0, &v, &sw, &st0](const auto& a, const auto& w, const auto& sv) {
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

  stan::test::expect_ad(f_a_w_sv(y, t0, v, sw, st0), a, w, sv);
}

TEST(mathMixScalFun_a_w_sw, wiener_full_lcdf) {
  auto f_a_w_sw = [](const auto& y, const auto& t0, const auto& v,
                     const auto& sv, const auto& st0) {
    return
        [&y, &t0, &v, &sv, &st0](const auto& a, const auto& w, const auto& sw) {
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

  stan::test::expect_ad(f_a_w_sw(y, t0, v, sv, st0), a, w, sw);
}

TEST(mathMixScalFun_a_w_st0, wiener_full_lcdf) {
  auto f_a_w_st0 = [](const auto& y, const auto& t0, const auto& v,
                      const auto& sv, const auto& sw) {
    return
        [&y, &t0, &v, &sv, &sw](const auto& a, const auto& w, const auto& st0) {
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

  stan::test::expect_ad(f_a_w_st0(y, t0, v, sv, sw), a, w, st0);
}

TEST(mathMixScalFun_a_v_sv, wiener_full_lcdf) {
  auto f_a_v_sv = [](const auto& y, const auto& t0, const auto& w,
                     const auto& sw, const auto& st0) {
    return
        [&y, &t0, &w, &sw, &st0](const auto& a, const auto& v, const auto& sv) {
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

  stan::test::expect_ad(f_a_v_sv(y, t0, w, sw, st0), a, v, sv);
}
