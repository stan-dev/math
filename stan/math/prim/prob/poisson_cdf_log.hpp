#ifndef STAN_MATH_PRIM_PROB_POISSON_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/poisson_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>poisson_lcdf</code>
 */
template <typename T_n, typename T_rate>
return_type_t<T_rate> poisson_cdf_log(const T_n& n, const T_rate& lambda) {
  return poisson_lcdf<T_n, T_rate>(n, lambda);
}

}  // namespace math
}  // namespace stan
#endif
