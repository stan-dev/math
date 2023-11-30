#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_y_a_w, wiener_full_lcdf) {
  auto f_y_a_w = [](const auto& t0, const auto& v, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&t0, &v, &sv, &sw, &st0](const auto& y, const auto& a, const auto& w) {
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

  stan::test::expect_ad(f_y_a_w(t0, v, sv, sw, st0), y, a, w);
}

TEST(mathMixScalFun_y_a_v, wiener_full_lcdf) {
  auto f_y_a_v = [](const auto& t0, const auto& w, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&t0, &w, &sv, &sw, &st0](const auto& y, const auto& a, const auto& v) {
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

  stan::test::expect_ad(f_y_a_v(t0, w, sv, sw, st0), y, a, v);
}

TEST(mathMixScalFun_y_a_sv, wiener_full_lcdf) {
  auto f_y_a_sv = [](const auto& t0, const auto& w, const auto& v,
                     const auto& sw, const auto& st0) {
    return
        [&t0, &w, &v, &sw, &st0](const auto& y, const auto& a, const auto& sv) {
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

  stan::test::expect_ad(f_y_a_sv(t0, w, v, sw, st0), y, a, sv);
}

TEST(mathMixScalFun_y_a_sw, wiener_full_lcdf) {
  auto f_y_a_sw = [](const auto& t0, const auto& w, const auto& v,
                     const auto& sv, const auto& st0) {
    return
        [&t0, &w, &v, &sv, &st0](const auto& y, const auto& a, const auto& sw) {
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

  stan::test::expect_ad(f_y_a_sw(t0, w, v, sv, st0), y, a, sw);
}

TEST(mathMixScalFun_y_a_st0, wiener_full_lcdf) {
  auto f_y_a_st0 = [](const auto& t0, const auto& w, const auto& v,
                      const auto& sv, const auto& sw) {
    return
        [&t0, &w, &v, &sv, &sw](const auto& y, const auto& a, const auto& st0) {
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

  stan::test::expect_ad(f_y_a_st0(t0, w, v, sv, sw), y, a, st0);
}
