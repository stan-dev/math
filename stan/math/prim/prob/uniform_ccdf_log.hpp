#ifndef STAN_MATH_PRIM_PROB_UNIFORM_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_UNIFORM_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/uniform_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>uniform_lccdf</code>
 */
template <typename T_y, typename T_low, typename T_high>
return_type_t<T_y, T_low, T_high> uniform_ccdf_log(const T_y& y,
                                                   const T_low& alpha,
                                                   const T_high& beta) {
  return uniform_lccdf<T_y, T_low, T_high>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
