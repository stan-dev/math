#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/lognormal_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>lognormal_lpdf</code>
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> lognormal_log(const T_y& y, const T_loc& mu,
                                                 const T_scale& sigma) {
  return lognormal_lpdf<propto, T_y, T_loc, T_scale>(y, mu, sigma);
}

/** \ingroup prob_dists
 * @deprecated use <code>lognormal_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> lognormal_log(const T_y& y,
                                                        const T_loc& mu,
                                                        const T_scale& sigma) {
  return lognormal_lpdf<T_y, T_loc, T_scale>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
