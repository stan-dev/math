#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_a_v_sw) {
  auto f_a_v_sw = [](const auto& y, const auto& t0, const auto& w,
                     const auto& sv, const auto& st0) {
    return
        [&y, &t0, &w, &sv, &st0](const auto& a, const auto& v, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_v_sw(y, t0, w, sv, st0), a, v, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_v_st0) {
  auto f_a_v_st0 = [](const auto& y, const auto& t0, const auto& w,
                      const auto& sv, const auto& sw) {
    return
        [&y, &t0, &w, &sv, &sw](const auto& a, const auto& v, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_v_st0(y, t0, w, sv, sw), a, v, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_sv_sw) {
  auto f_a_sv_sw = [](const auto& y, const auto& t0, const auto& w,
                      const auto& v, const auto& st0) {
    return
        [&y, &t0, &w, &v, &st0](const auto& a, const auto& sv, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_sv_sw(y, t0, w, v, st0), a, sv, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_sv_st0) {
  auto f_a_sv_st0 = [](const auto& y, const auto& t0, const auto& w,
                       const auto& v, const auto& sw) {
    return
        [&y, &t0, &w, &v, &sw](const auto& a, const auto& sv, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_sv_st0(y, t0, w, v, sw), a, sv, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_a_sw_st0) {
  auto f_a_sw_st0 = [](const auto& y, const auto& t0, const auto& w,
                       const auto& v, const auto& sv) {
    return
        [&y, &t0, &w, &v, &sv](const auto& a, const auto& sw, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_sw_st0(y, t0, w, v, sv), a, sw, st0);
}
