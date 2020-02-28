#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/exponential_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>exponential_lpdf</code>
 */
template <bool propto, typename T_y, typename T_inv_scale>
return_type_t<T_y, T_inv_scale> exponential_log(const T_y& y,
                                                const T_inv_scale& beta) {
  return exponential_lpdf<propto, T_y, T_inv_scale>(y, beta);
}

/** \ingroup prob_dists
 * @deprecated use <code>exponential_lpdf</code>
 */
template <typename T_y, typename T_inv_scale>
inline return_type_t<T_y, T_inv_scale> exponential_log(
    const T_y& y, const T_inv_scale& beta) {
  return exponential_lpdf<T_y, T_inv_scale>(y, beta);
}

}  // namespace math
}  // namespace stan
#endif
