#ifndef STAN_MATH_PRIM_PROB_WEIBULL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_WEIBULL_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/weibull_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>weibull_lccdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> weibull_ccdf_log(const T_y& y,
                                                      const T_shape& alpha,
                                                      const T_scale& sigma) {
  return weibull_lccdf<T_y, T_shape, T_scale>(y, alpha, sigma);
}

}  // namespace math
}  // namespace stan
#endif
