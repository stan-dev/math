#ifndef STAN_MATH_PRIM_PROB_VON_MISES_LOG_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_LOG_HPP

#include <stanh/prim/meta/return_type.hpp>
#include <stanh/prim/prob/von_mises_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>von_mises_lpdf</code>
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
typename return_type<T_y, T_loc, T_scale>::type von_mises_log(
    T_y const& y, T_loc const& mu, T_scale const& kappa) {
  return von_mises_lpdf<propto, T_y, T_loc, T_scale>(y, mu, kappa);
}

/**
 * @deprecated use <code>von_mises_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline typename return_type<T_y, T_loc, T_scale>::type von_mises_log(
    T_y const& y, T_loc const& mu, T_scale const& kappa) {
  return von_mises_lpdf<T_y, T_loc, T_scale>(y, mu, kappa);
}

}  // namespace math
}  // namespace stan
#endif
