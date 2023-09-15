#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv> wiener5_lpdf(
    const T_y& y, const T_a& a, const T_t0& t0, const T_w& w, const T_v& v,
    const T_sv& sv, const double& precision_derivatives);

TEST(mathMixScalFun, wiener5_lpdf) {
  using stan::math::fvar;
  using stan::math::var;
  double y = 0;
  double a = 0;
  double t0 = 0;
  double w = 0;
  double v = 0;
  double sv = 0;
  double sw = 0;
  double st0 = 0;
  stan::math::wiener5_lpdf(y, a, t0, w, v, sv, st0);

}
