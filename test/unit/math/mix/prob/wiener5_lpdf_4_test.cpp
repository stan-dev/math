#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener5MixArgs, wiener_lpdf_a_v_sv) {
  auto f_a_v_sv = [](const auto& y, const auto& t0, const auto& w) {
    return [&y, &t0, &w](const auto& a, const auto& v, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_v_sv(y, t0, w), a, v, sv);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_t0_w_v) {
  auto f_t0_w_v = [](const auto& y, const auto& a, const auto& sv) {
    return [&y, &a, &sv](const auto& t0, const auto& w, const auto& v) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_t0_w_v(y, a, sv), t0, w, v);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_t0_w_sv) {
  auto f_t0_w_sv = [](const auto& y, const auto& a, const auto& v) {
    return [&y, &a, &v](const auto& t0, const auto& w, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_t0_w_sv(y, a, v), t0, w, sv);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_t0_v_sv) {
  auto f_t0_v_sv = [](const auto& y, const auto& a, const auto& w) {
    return [&y, &a, &w](const auto& t0, const auto& v, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_t0_v_sv(y, a, w), t0, v, sv);
}

TEST_F(Wiener5MixArgs, wiener_lpdf_w_v_sv) {
  auto f_w_v_sv = [](const auto& y, const auto& a, const auto& t0) {
    return [&y, &a, &t0](const auto& w, const auto& v, const auto& sv) {
      return stan::math::wiener_lpdf(y, a, t0, w, v, sv, 1e-4);
    };
  };
  stan::test::expect_ad(f_w_v_sv(y, a, t0), w, v, sv);
}
