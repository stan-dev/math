#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_w_v_sv) {
  auto f_w_v_sv = [](const auto& y, const auto& a, const auto& t0,
                     const auto& sw, const auto& st0) {
    return
        [&y, &a, &t0, &sw, &st0](const auto& w, const auto& v, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_w_v_sv(y, a, t0, sw, st0), w, v, sv);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_w_v_sw) {
  auto f_w_v_sw = [](const auto& y, const auto& a, const auto& t0,
                     const auto& sv, const auto& st0) {
    return
        [&y, &a, &t0, &sv, &st0](const auto& w, const auto& v, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_w_v_sw(y, a, t0, sv, st0), w, v, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_w_v_st0) {
  auto f_w_v_st0 = [](const auto& y, const auto& a, const auto& t0,
                      const auto& sv, const auto& sw) {
    return
        [&y, &a, &t0, &sv, &sw](const auto& w, const auto& v, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_w_v_st0(y, a, t0, sv, sw), w, v, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_w_sv_sw) {
  auto f_w_sv_sw = [](const auto& y, const auto& a, const auto& t0,
                      const auto& v, const auto& st0) {
    return
        [&y, &a, &t0, &v, &st0](const auto& w, const auto& sv, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_w_sv_sw(y, a, t0, v, st0), w, sv, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_w_sv_st0) {
  auto f_w_sv_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& v, const auto& sw) {
    return
        [&y, &a, &t0, &v, &sw](const auto& w, const auto& sv, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_w_sv_st0(y, a, t0, v, sw), w, sv, st0);
}
