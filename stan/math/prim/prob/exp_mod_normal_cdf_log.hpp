#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/exp_mod_normal_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>exp_mod_normal_lcdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_cdf_log(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  return exp_mod_normal_lcdf<T_y, T_loc, T_scale, T_inv_scale>(y, mu, sigma,
                                                               lambda);
}

}  // namespace math
}  // namespace stan
#endif
