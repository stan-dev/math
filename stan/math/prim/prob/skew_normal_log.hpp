#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/skew_normal_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>skew_normal_lpdf</code>
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_log(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  return skew_normal_lpdf<propto, T_y, T_loc, T_scale, T_shape>(y, mu, sigma,
                                                                alpha);
}

/** \ingroup prob_dists
 * @deprecated use <code>skew_normal_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_log(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  return skew_normal_lpdf<T_y, T_loc, T_scale, T_shape>(y, mu, sigma, alpha);
}

}  // namespace math
}  // namespace stan
#endif
