#ifndef STAN_MATH_PRIM_PROB_MULTI_NORMAL_PREC_LOG_HPP
#define STAN_MATH_PRIM_PROB_MULTI_NORMAL_PREC_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/multi_normal_prec_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>multi_normal_prec_lpdf</code>
 */
template <bool propto, typename T_y, typename T_loc, typename T_covar>
return_type_t<T_y, T_loc, T_covar> multi_normal_prec_log(const T_y& y,
                                                         const T_loc& mu,
                                                         const T_covar& Sigma) {
  return multi_normal_prec_lpdf<propto, T_y, T_loc, T_covar>(y, mu, Sigma);
}

/** \ingroup multivar_dists
 * @deprecated use <code>multi_normal_prec_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_covar>
inline return_type_t<T_y, T_loc, T_covar> multi_normal_prec_log(
    const T_y& y, const T_loc& mu, const T_covar& Sigma) {
  return multi_normal_prec_lpdf<T_y, T_loc, T_covar>(y, mu, Sigma);
}

}  // namespace math
}  // namespace stan
#endif
