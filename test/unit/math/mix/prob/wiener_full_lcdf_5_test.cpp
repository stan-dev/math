#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/util.hpp>

TEST_F(Wiener7MixArgs, wiener_lcdf_a_t0_w) {
  auto f_a_t0_w = [](const auto& y, const auto& v, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&y, &v, &sv, &sw, &st0](const auto& a, const auto& t0, const auto& w) {
          return stan::math::wiener_lcdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_t0_w(y, v, sv, sw, st0), a, t0, w);
}

TEST_F(Wiener7MixArgs, wiener_lcdf_a_t0_v) {
  auto f_a_t0_v = [](const auto& y, const auto& w, const auto& sv,
                     const auto& sw, const auto& st0) {
    return
        [&y, &w, &sv, &sw, &st0](const auto& a, const auto& t0, const auto& v) {
          return stan::math::wiener_lcdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_t0_v(y, w, sv, sw, st0), a, t0, v);
}

TEST_F(Wiener7MixArgs, wiener_lcdf_a_t0_sv) {
  auto f_a_t0_sv = [](const auto& y, const auto& w, const auto& v,
                      const auto& sw, const auto& st0) {
    return
        [&y, &w, &v, &sw, &st0](const auto& a, const auto& t0, const auto& sv) {
          return stan::math::wiener_lcdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_t0_sv(y, w, v, sw, st0), a, t0, sv);
}

TEST_F(Wiener7MixArgs, wiener_lcdf_a_t0_sw) {
  auto f_a_t0_sw = [](const auto& y, const auto& w, const auto& v,
                      const auto& sv, const auto& st0) {
    return
        [&y, &w, &v, &sv, &st0](const auto& a, const auto& t0, const auto& sw) {
          return stan::math::wiener_lcdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_t0_sw(y, w, v, sv, st0), a, t0, sw);
}

TEST_F(Wiener7MixArgs, wiener_lcdf_a_t0_st0) {
  auto f_a_t0_st0 = [](const auto& y, const auto& w, const auto& v,
                       const auto& sv, const auto& sw) {
    return
        [&y, &w, &v, &sv, &sw](const auto& a, const auto& t0, const auto& st0) {
          return stan::math::wiener_lcdf(y, a, t0, w, v, sv, sw, st0);
        };
  };
  stan::test::expect_ad(f_a_t0_st0(y, w, v, sv, sw), a, t0, st0);
}
