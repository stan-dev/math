#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun_a_t0_w, wiener5_lpdf) {
  auto f_a_t0_w = [](const auto& y, const auto& v, const auto& sv) {
    return [&y, &v, &sv](const auto& a, const auto& t0, const auto& w) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.5;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_t0_w(y, v, sv), a, t0, w);
}

TEST(mathMixScalFun_a_t0_v, wiener5_lpdf) {
  auto f_a_t0_v = [](const auto& y, const auto& w, const auto& sv) {
    return [&y, &w, &sv](const auto& a, const auto& t0, const auto& v) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.5;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_t0_v(y, w, sv), a, t0, v);
}

TEST(mathMixScalFun_a_t0_sv, wiener5_lpdf) {
  auto f_a_t0_sv = [](const auto& y, const auto& w, const auto& v) {
    return [&y, &w, &v](const auto& a, const auto& t0, const auto& sv) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.5;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_t0_sv(y, w, v), a, t0, sv);
}

TEST(mathMixScalFun_a_w_v, wiener5_lpdf) {
  auto f_a_w_v = [](const auto& y, const auto& t0, const auto& sv) {
    return [&y, &t0, &sv](const auto& a, const auto& w, const auto& v) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.5;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_w_v(y, t0, sv), a, w, v);
}

TEST(mathMixScalFun_a_w_sv, wiener5_lpdf) {
  auto f_a_w_sv = [](const auto& y, const auto& t0, const auto& v) {
    return [&y, &t0, &v](const auto& a, const auto& w, const auto& sv) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.5;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_w_sv(y, t0, v), a, w, sv);
}
