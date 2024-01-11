#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_a_w_v) {
  auto f_a_w_v = [](const auto& y, const auto& t0, const auto& sv,
                    const auto& sw, const auto& st0) {
    return
        [&y, &t0, &sv, &sw, &st0](const auto& a, const auto& w, const auto& v) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_w_v(y, t0, sv, sw, st0), a, w, v);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_w_sv) {
  auto f_a_w_sv = [](const auto& y, const auto& t0, const auto& v,
                     const auto& sw, const auto& st0) {
    return
        [&y, &t0, &v, &sw, &st0](const auto& a, const auto& w, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_w_sv(y, t0, v, sw, st0), a, w, sv);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_w_sw) {
  auto f_a_w_sw = [](const auto& y, const auto& t0, const auto& v,
                     const auto& sv, const auto& st0) {
    return
        [&y, &t0, &v, &sv, &st0](const auto& a, const auto& w, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_w_sw(y, t0, v, sv, st0), a, w, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_w_st0) {
  auto f_a_w_st0 = [](const auto& y, const auto& t0, const auto& v,
                      const auto& sv, const auto& sw) {
    return
        [&y, &t0, &v, &sv, &sw](const auto& a, const auto& w, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_w_st0(y, t0, v, sv, sw), a, w, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_v_sv) {
  auto f_a_v_sv = [](const auto& y, const auto& t0, const auto& w,
                     const auto& sw, const auto& st0) {
    return
        [&y, &t0, &w, &sw, &st0](const auto& a, const auto& v, const auto& sv) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_v_sv(y, t0, w, sw, st0), a, v, sv);
}
