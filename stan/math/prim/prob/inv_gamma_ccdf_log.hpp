#ifndef STAN_MATH_PRIM_PROB_INV_GAMMA_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_INV_GAMMA_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/inv_gamma_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>inv_gamma_lccdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> inv_gamma_ccdf_log(const T_y& y,
                                                        const T_shape& alpha,
                                                        const T_scale& beta) {
  return inv_gamma_lccdf<T_y, T_shape, T_scale>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
