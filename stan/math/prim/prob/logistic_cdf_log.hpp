#ifndef STAN_MATH_PRIM_PROB_LOGISTIC_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_LOGISTIC_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/logistic_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>logistic_lcdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> logistic_cdf_log(const T_y& y,
                                                    const T_loc& mu,
                                                    const T_scale& sigma) {
  return logistic_lcdf<T_y, T_loc, T_scale>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
