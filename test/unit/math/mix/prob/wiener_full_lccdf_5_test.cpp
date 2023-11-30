#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_a_t0_w, wiener_full_lccdf) {
  auto f_a_t0_w = [](const auto& y, const auto& v, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&y, &v, &sv, &sw, &st0](const auto& a, const auto& t0, const auto& w) {
          return stan::math::wiener_full_lccdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_a_t0_w(y, v, sv, sw, st0), a, t0, w);
}

TEST(mathMixScalFun_a_t0_v, wiener_full_lccdf) {
  auto f_a_t0_v = [](const auto& y, const auto& w, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&y, &w, &sv, &sw, &st0](const auto& a, const auto& t0, const auto& v) {
          return stan::math::wiener_full_lccdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_a_t0_v(y, w, sv, sw, st0), a, t0, v);
}

TEST(mathMixScalFun_a_t0_sv, wiener_full_lccdf) {
  auto f_a_t0_sv = [](const auto& y, const auto& w, const auto& v,
                      const auto& sw, const auto& st0) {
    return
        [&y, &w, &v, &sw, &st0](const auto& a, const auto& t0, const auto& sv) {
          return stan::math::wiener_full_lccdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_a_t0_sv(y, w, v, sw, st0), a, t0, sv);
}

TEST(mathMixScalFun_a_t0_sw, wiener_full_lccdf) {
  auto f_a_t0_sw = [](const auto& y, const auto& w, const auto& v,
                      const auto& sv, const auto& st0) {
    return
        [&y, &w, &v, &sv, &st0](const auto& a, const auto& t0, const auto& sw) {
          return stan::math::wiener_full_lccdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_a_t0_sw(y, w, v, sv, st0), a, t0, sw);
}

TEST(mathMixScalFun_a_t0_st0, wiener_full_lccdf) {
  auto f_a_t0_st0 = [](const auto& y, const auto& w, const auto& v,
                       const auto& sv, const auto& sw) {
    return
        [&y, &w, &v, &sv, &sw](const auto& a, const auto& t0, const auto& st0) {
          return stan::math::wiener_full_lccdf(y, a, t0, w, v, sv, sw, st0);
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

  stan::test::expect_ad(f_a_t0_st0(y, w, v, sv, sw), a, t0, st0);
}
