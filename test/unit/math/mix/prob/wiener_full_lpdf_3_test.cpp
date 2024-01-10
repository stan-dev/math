#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_y_w_v) {
  auto f_y_w_v = [](const auto& a, const auto& t0, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&a, &t0, &sv, &sw, &st0](const auto& y, const auto& w, const auto& v) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_w_v(a, t0, sv, sw, st0), y, w, v);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_w_sv) {
  auto f_y_w_sv = [](const auto& a, const auto& t0, const auto& v,
                     const auto& sw, const auto& st0) {
    return
        [&a, &t0, &v, &sw, &st0](const auto& y, const auto& w, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_w_sv(a, t0, v, sw, st0), y, w, sv);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_w_sw) {
  auto f_y_w_sw = [](const auto& a, const auto& t0, const auto& v,
                     const auto& sv, const auto& st0) {
    return
        [&a, &t0, &v, &sv, &st0](const auto& y, const auto& w, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_w_sw(a, t0, v, sv, st0), y, w, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_w_st0) {
  auto f_y_w_st0 = [](const auto& a, const auto& t0, const auto& v,
                      const auto& sv, const auto& sw) {
    return
        [&a, &t0, &v, &sv, &sw](const auto& y, const auto& w, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_w_st0(a, t0, v, sv, sw), y, w, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_v_sv) {
  auto f_y_v_sv = [](const auto& a, const auto& t0, const auto& w,
                     const auto& sw, const auto& st0) {
    return
        [&a, &t0, &w, &sw, &st0](const auto& y, const auto& v, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_v_sv(a, t0, w, sw, st0), y, v, sv);
}
