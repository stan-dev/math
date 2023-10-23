#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_y_a_t0, wiener5_lpdf) {
  auto f_y_a_t0 = [](const auto& w, const auto& v, const auto& sv) {
    return [&w, &v, &sv](const auto& y, const auto& a, const auto& t0) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_y_a_t0(w, v, sv), y, a, t0);
}

TEST(mathMixScalFun_y_a_w, wiener5_lpdf) {
  auto f_y_a_w = [](const auto& t0, const auto& v, const auto& sv) {
    return [&t0, &v, &sv](const auto& y, const auto& a, const auto& w) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_y_a_w(t0, v, sv), y, a, w);
}

TEST(mathMixScalFun_y_a_v, wiener5_lpdf) {
  auto f_y_a_v = [](const auto& t0, const auto& w, const auto& sv) {
    return [&t0, &w, &sv](const auto& y, const auto& a, const auto& v) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_y_a_v(t0, w, sv), y, a, v);
}

TEST(mathMixScalFun_y_a_sv, wiener5_lpdf) {
  auto f_y_a_sv = [](const auto& t0, const auto& w, const auto& v) {
    return [&t0, &w, &v](const auto& y, const auto& a, const auto& sv) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_y_a_sv(t0, w, v), y, a, sv);
}

TEST(mathMixScalFun_y_t0_w, wiener5_lpdf) {
  auto f_y_t0_w = [](const auto& a, const auto& v, const auto& sv) {
    return [&a, &v, &sv](const auto& y, const auto& t0, const auto& w) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_y_t0_w(a, v, sv), y, t0, w);
}
