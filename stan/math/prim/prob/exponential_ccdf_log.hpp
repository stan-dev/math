#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/exponential_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>exponential_lccdf</code>
 */
template <typename T_y, typename T_inv_scale>
return_type_t<T_y, T_inv_scale> exponential_ccdf_log(const T_y& y,
                                                     const T_inv_scale& beta) {
  return exponential_lccdf<T_y, T_inv_scale>(y, beta);
}

}  // namespace math
}  // namespace stan
#endif
