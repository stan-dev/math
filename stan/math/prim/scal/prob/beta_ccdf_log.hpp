#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/beta_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>beta_lccdf</code>
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
inline auto beta_ccdf_log(T_y&& y, T_scale_succ&& alpha, T_scale_fail&& beta) {
  return beta_lccdf(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
