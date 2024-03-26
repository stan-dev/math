#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener5MixArgs, wiener_lpdf_a_t0_w) {
  auto f_a_t0_w = [](const auto& y, const auto& v, const auto& sv) {
    return [&y, &v, &sv](const auto& a, const auto& t0, const auto& w) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_t0_w(y, v, sv), a, t0, w);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_a_t0_v) {
  auto f_a_t0_v = [](const auto& y, const auto& w, const auto& sv) {
    return [&y, &w, &sv](const auto& a, const auto& t0, const auto& v) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_t0_v(y, w, sv), a, t0, v);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_a_t0_sv) {
  auto f_a_t0_sv = [](const auto& y, const auto& w, const auto& v) {
    return [&y, &w, &v](const auto& a, const auto& t0, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_t0_sv(y, w, v), a, t0, sv);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_a_w_v) {
  auto f_a_w_v = [](const auto& y, const auto& t0, const auto& sv) {
    return [&y, &t0, &sv](const auto& a, const auto& w, const auto& v) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_w_v(y, t0, sv), a, w, v);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_a_w_sv) {
  auto f_a_w_sv = [](const auto& y, const auto& t0, const auto& v) {
    return [&y, &t0, &v](const auto& a, const auto& w, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_w_sv(y, t0, v), a, w, sv);
}
