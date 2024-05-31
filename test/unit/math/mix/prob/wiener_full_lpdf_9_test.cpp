#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_v_sw) {
  auto f_t0_v_sw = [](const auto& y, const auto& a, const auto& w,
                      const auto& sv, const auto& st0) {
    return
        [&y, &a, &w, &sv, &st0](const auto& t0, const auto& v, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_v_sw(y, a, w, sv, st0), t0, v, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_v_st0) {
  auto f_t0_v_st0 = [](const auto& y, const auto& a, const auto& w,
                       const auto& sv, const auto& sw) {
    return
        [&y, &a, &w, &sv, &sw](const auto& t0, const auto& v, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_v_st0(y, a, w, sv, sw), t0, v, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_sv_sw) {
  auto f_t0_sv_sw = [](const auto& y, const auto& a, const auto& w,
                       const auto& v, const auto& st0) {
    return
        [&y, &a, &w, &v, &st0](const auto& t0, const auto& sv, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_sv_sw(y, a, w, v, st0), t0, sv, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_sv_st0) {
  auto f_t0_sv_st0 = [](const auto& y, const auto& a, const auto& w,
                        const auto& v, const auto& sw) {
    return
        [&y, &a, &w, &v, &sw](const auto& t0, const auto& sv, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_sv_st0(y, a, w, v, sw), t0, sv, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_t0_sw_st0) {
  auto f_t0_sw_st0 = [](const auto& y, const auto& a, const auto& w,
                        const auto& v, const auto& sv) {
    return
        [&y, &a, &w, &v, &sv](const auto& t0, const auto& sw, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_t0_sw_st0(y, a, w, v, sv), t0, sw, st0);
}
