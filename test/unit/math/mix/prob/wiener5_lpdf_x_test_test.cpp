#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

/*TEST(mathMixScalFun_a_v_sv, wiener5_lpdf) {
  auto f_a_v_sv = [](const auto& y, const auto& t0, const auto& w) {
    return [&y, &t0, &w](const auto& a, const auto& v, const auto& sv) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_v_sv(y, t0, w), a, v, sv);
}*/

/*TEST(mathMixScalFun_a_v_sv, wiener5_lpdf) {
  auto f_a_v_sv = [](const auto& y, const auto& t0, const auto& w, const auto& sv) {
    return [&y, &t0, &w, &sv](const auto& a, const auto& v) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_v_sv(y, t0, w, sv), a, v);
}*/

TEST(mathMixScalFun_a_v_sv, wiener5_lpdf) {
  auto f_a_v_sv = [](const auto& y, const auto& t0, const auto& w, const auto& v, const auto& sv) {
    return [&y, &t0, &w, &v, &sv](const auto& a) {
      return stan::math::wiener5_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };

  double y = 1.0;
  double a = 2.0;
  double t0 = 0.2;
  double w = 0.5;
  double v = 2.0;
  double sv = 0.2;

  stan::test::expect_ad(f_a_v_sv(y, t0, w, v, sv), a);
}


