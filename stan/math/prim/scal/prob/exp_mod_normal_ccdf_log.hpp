#ifndef STAN_MATH_PRIM_SCAL_PROB_EXP_MOD_NORMAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXP_MOD_NORMAL_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/exp_mod_normal_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>exp_mod_normal_lccdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
inline auto exp_mod_normal_ccdf_log(const T_y& y, const T_loc& mu,
                                    const T_scale& sigma,
                                    const T_inv_scale& lambda) {
  return exp_mod_normal_lccdf(y, mu, sigma, lambda);
}

}  // namespace math
}  // namespace stan
#endif
