#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lpdf_y_v_sw) {
  auto f_y_v_sw = [](const auto& a, const auto& t0, const auto& w,
                     const auto& sv, const auto& st0) {
    return
        [&a, &t0, &w, &sv, &st0](const auto& y, const auto& v, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_v_sw(a, t0, w, sv, st0), y, v, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_v_st0) {
  auto f_y_v_st0 = [](const auto& a, const auto& t0, const auto& w,
                      const auto& sv, const auto& sw) {
    return
        [&a, &t0, &w, &sv, &sw](const auto& y, const auto& v, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_v_st0(a, t0, w, sv, sw), y, v, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_sv_sw) {
  auto f_y_sv_sw = [](const auto& a, const auto& t0, const auto& w,
                      const auto& v, const auto& st0) {
    return
        [&a, &t0, &w, &v, &st0](const auto& y, const auto& sv, const auto& sw) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_sv_sw(a, t0, w, v, st0), y, sv, sw);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_sv_st0) {
  auto f_y_sv_st0 = [](const auto& a, const auto& t0, const auto& w,
                       const auto& v, const auto& sw) {
    return
        [&a, &t0, &w, &v, &sw](const auto& y, const auto& sv, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_sv_st0(a, t0, w, v, sw), y, sv, st0);
}

TEST_F(Wiener7MixArgs, wiener_lpdf_y_sw_st0) {
  auto f_y_sw_st0 = [](const auto& a, const auto& t0, const auto& w,
                       const auto& v, const auto& sv) {
    return
        [&a, &t0, &w, &v, &sv](const auto& y, const auto& sw, const auto& st0) {
          return stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_y_sw_st0(a, t0, w, v, sv), y, sw, st0);
}
