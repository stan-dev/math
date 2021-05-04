#ifndef STAN_MATH_PRIM_PROB_DISCRETE_RANGE_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_DISCRETE_RANGE_CDF_LOG_HPP

#include <stan/math/prim/prob/discrete_range_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>discrete_range_lcdf</code>
 */
template <typename T_y, typename T_lower, typename T_upper>
double discrete_range_cdf_log(const T_y& y, const T_lower& lower,
                              const T_upper& upper) {
  return discrete_range_lcdf(y, lower, upper);
}

}  // namespace math
}  // namespace stan
#endif
