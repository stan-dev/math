#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener5MixArgs, wiener_lpdf_y_t0_v) {
  auto f_y_t0_v = [](const auto& a, const auto& w, const auto& sv) {
    return [&a, &w, &sv](const auto& y, const auto& t0, const auto& v) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_t0_v(a, w, sv), y, t0, v);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_y_t0_sv) {
  auto f_y_t0_sv = [](const auto& a, const auto& w, const auto& v) {
    return [&a, &w, &v](const auto& y, const auto& t0, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_t0_sv(a, w, v), y, t0, sv);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_y_w_v) {
  auto f_y_w_v = [](const auto& a, const auto& t0, const auto& sv) {
    return [&a, &t0, &sv](const auto& y, const auto& w, const auto& v) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_w_v(a, t0, sv), y, w, v);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_y_w_sv) {
  auto f_y_w_sv = [](const auto& a, const auto& t0, const auto& v) {
    return [&a, &t0, &v](const auto& y, const auto& w, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_w_sv(a, t0, v), y, w, sv);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_y_v_sv) {
  auto f_y_v_sv = [](const auto& a, const auto& t0, const auto& w) {
    return [&a, &t0, &w](const auto& y, const auto& v, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_v_sv(a, t0, w), y, v, sv);
}
