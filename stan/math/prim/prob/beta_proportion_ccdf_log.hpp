#ifndef STAN_MATH_PRIM_PROB_BETA_PROPORTION_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_BETA_PROPORTION_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/beta_proportion_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>beta_proportion_lccdf</code>
 */
template <typename T_y, typename T_loc, typename T_prec>
return_type_t<T_y, T_loc, T_prec> beta_proportion_ccdf_log(
    const T_y& y, const T_loc& mu, const T_prec& kappa) {
  return beta_proportion_lccdf<T_y, T_loc, T_prec>(y, mu, kappa);
}

}  // namespace math
}  // namespace stan
#endif
