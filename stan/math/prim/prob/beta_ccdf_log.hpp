#ifndef STAN_MATH_PRIM_PROB_BETA_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_BETA_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/beta_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>beta_lccdf</code>
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
return_type_t<T_y, T_scale_succ, T_scale_fail> beta_ccdf_log(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta) {
  return beta_lccdf<T_y, T_scale_succ, T_scale_fail>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
