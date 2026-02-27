#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_y_a_t0) {
  auto f_y_a_t0 = [](const auto& w, const auto& v) {
    return [&w, &v](const auto& y, const auto& a, const auto& t0) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_a_t0(w, v), y, a, t0);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_y_a_w) {
  auto f_y_a_w = [](const auto& t0, const auto& v) {
    return [&t0, &v](const auto& y, const auto& a, const auto& w) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_a_w(t0, v), y, a, w);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_y_a_v) {
  auto f_y_a_v = [](const auto& t0, const auto& w) {
    return [&t0, &w](const auto& y, const auto& a, const auto& v) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_a_v(t0, w), y, a, v);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_y_t0_w) {
  auto f_y_t0_w = [](const auto& a, const auto& v) {
    return [&a, &v](const auto& y, const auto& t0, const auto& w) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_t0_w(a, v), y, t0, w);
}

TEST_F(Wiener4MixArgs, wiener_lccdf_unnorm_y_t0_v) {
  auto f_y_t0_v = [](const auto& a, const auto& w) {
    return [&a, &w](const auto& y, const auto& t0, const auto& v) {
      return stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, 1e-4);
    };
  };
  stan::test::expect_ad(f_y_t0_v(a, w), y, t0, v);
}
