#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_y_t0_w) {
  auto f_y_t0_w = [](const auto& a, const auto& v, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&a, &v, &sv, &sw, &st0](const auto& y, const auto& t0, const auto& w) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_t0_w(a, v, sv, sw, st0), y, t0, w);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_t0_v) {
  auto f_y_t0_v = [](const auto& a, const auto& w, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&a, &w, &sv, &sw, &st0](const auto& y, const auto& t0, const auto& v) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_t0_v(a, w, sv, sw, st0), y, t0, v);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_t0_sv) {
  auto f_y_t0_sv = [](const auto& a, const auto& w, const auto& v,
                      const auto& sw, const auto& st0) {
    return
        [&a, &w, &v, &sw, &st0](const auto& y, const auto& t0, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_t0_sv(a, w, v, sw, st0), y, t0, sv);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_t0_sw) {
  auto f_y_t0_sw = [](const auto& a, const auto& w, const auto& v,
                      const auto& sv, const auto& st0) {
    return
        [&a, &w, &v, &sv, &st0](const auto& y, const auto& t0, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_t0_sw(a, w, v, sv, st0), y, t0, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_t0_st0) {
  auto f_y_t0_st0 = [](const auto& a, const auto& w, const auto& v,
                       const auto& sv, const auto& sw) {
    return
        [&a, &w, &v, &sv, &sw](const auto& y, const auto& t0, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_t0_st0(a, w, v, sv, sw), y, t0, st0);
}
