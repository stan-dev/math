#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_y_w_v) {
  auto f_y_w_v = [](const auto& a, const auto& t0) {
    return [&a, &t0](const auto& y, const auto& w, const auto& v) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_w_v(a, t0), y, w, v);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_a_t0_w) {
  auto f_a_t0_w = [](const auto& y, const auto& v) {
    return [&y, &v](const auto& a, const auto& t0, const auto& w) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_t0_w(y, v), a, t0, w);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_a_t0_v) {
  auto f_a_t0_v = [](const auto& y, const auto& w) {
    return [&y, &w](const auto& a, const auto& t0, const auto& v) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_t0_v(y, w), a, t0, v);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_a_w_v) {
  auto f_a_w_v = [](const auto& y, const auto& t0) {
    return [&y, &t0](const auto& a, const auto& w, const auto& v) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_a_w_v(y, t0), a, w, v);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_t0_w_v) {
  auto f_t0_w_v = [](const auto& y, const auto& a) {
    return [&y, &a](const auto& t0, const auto& w, const auto& v) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_t0_w_v(y, a), t0, w, v);
}
