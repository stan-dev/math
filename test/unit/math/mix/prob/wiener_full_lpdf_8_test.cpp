#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_w_v) {
  auto f_t0_w_v = [](const auto& y, const auto& a, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&y, &a, &sv, &sw, &st0](const auto& t0, const auto& w, const auto& v) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_w_v(y, a, sv, sw, st0), t0, w, v);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_w_sv) {
  auto f_t0_w_sv = [](const auto& y, const auto& a, const auto& v,
                      const auto& sw, const auto& st0) {
    return
        [&y, &a, &v, &sw, &st0](const auto& t0, const auto& w, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_w_sv(y, a, v, sw, st0), t0, w, sv);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_w_sw) {
  auto f_t0_w_sw = [](const auto& y, const auto& a, const auto& v,
                      const auto& sv, const auto& st0) {
    return
        [&y, &a, &v, &sv, &st0](const auto& t0, const auto& w, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_w_sw(y, a, v, sv, st0), t0, w, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_w_st0) {
  auto f_t0_w_st0 = [](const auto& y, const auto& a, const auto& v,
                       const auto& sv, const auto& sw) {
    return
        [&y, &a, &v, &sv, &sw](const auto& t0, const auto& w, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_w_st0(y, a, v, sv, sw), t0, w, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_v_sv) {
  auto f_t0_v_sv = [](const auto& y, const auto& a, const auto& w,
                      const auto& sw, const auto& st0) {
    return
        [&y, &a, &w, &sw, &st0](const auto& t0, const auto& v, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_v_sv(y, a, w, sw, st0), t0, v, sv);
}
