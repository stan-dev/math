#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_y_a_w) {
  auto f_y_a_w = [](const auto& t0, const auto& v, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&t0, &v, &sv, &sw, &st0](const auto& y, const auto& a, const auto& w) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_a_w(t0, v, sv, sw, st0), y, a, w);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_a_v) {
  auto f_y_a_v = [](const auto& t0, const auto& w, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&t0, &w, &sv, &sw, &st0](const auto& y, const auto& a, const auto& v) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_a_v(t0, w, sv, sw, st0), y, a, v);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_a_sv) {
  auto f_y_a_sv = [](const auto& t0, const auto& w, const auto& v,
                     const auto& sw, const auto& st0) {
    return
        [&t0, &w, &v, &sw, &st0](const auto& y, const auto& a, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_a_sv(t0, w, v, sw, st0), y, a, sv);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_a_sw) {
  auto f_y_a_sw = [](const auto& t0, const auto& w, const auto& v,
                     const auto& sv, const auto& st0) {
    return
        [&t0, &w, &v, &sv, &st0](const auto& y, const auto& a, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_a_sw(t0, w, v, sv, st0), y, a, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_a_st0) {
  auto f_y_a_st0 = [](const auto& t0, const auto& w, const auto& v,
                      const auto& sv, const auto& sw) {
    return
        [&t0, &w, &v, &sv, &sw](const auto& y, const auto& a, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_a_st0(t0, w, v, sv, sw), y, a, st0);
}
