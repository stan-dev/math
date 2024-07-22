#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_w_sw_st0) {
  auto f_w_sw_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& v, const auto& sv) {
    return
        [&y, &a, &t0, &v, &sv](const auto& w, const auto& sw, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_w_sw_st0(y, a, t0, v, sv), w, sw, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_v_sv_sw) {
  auto f_v_sv_sw = [](const auto& y, const auto& a, const auto& t0,
                      const auto& w, const auto& st0) {
    return
        [&y, &a, &t0, &w, &st0](const auto& v, const auto& sv, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_v_sv_sw(y, a, t0, w, st0), v, sv, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_v_sv_st0) {
  auto f_v_sv_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& w, const auto& sw) {
    return
        [&y, &a, &t0, &w, &sw](const auto& v, const auto& sv, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_v_sv_st0(y, a, t0, w, sw), v, sv, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_v_sw_st0) {
  auto f_v_sw_st0 = [](const auto& y, const auto& a, const auto& t0,
                       const auto& w, const auto& sv) {
    return
        [&y, &a, &t0, &w, &sv](const auto& v, const auto& sw, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_v_sw_st0(y, a, t0, w, sv), v, sw, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_sv_sw_st0) {
  auto f_sv_sw_st0 = [](const auto& y, const auto& a, const auto& t0,
                        const auto& w, const auto& v) {
    return
        [&y, &a, &t0, &w, &v](const auto& sv, const auto& sw, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_sv_sw_st0(y, a, t0, w, v), sv, sw, st0);
}
