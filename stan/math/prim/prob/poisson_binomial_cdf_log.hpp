#ifndef STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/poisson_binomial_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>poisson_binomial_lcdf</code>
 */
template <typename T_y, typename T_theta>
return_type_t<T_theta> poisson_binomial_cdf_log(const T_y& y,
                                                const T_theta& theta) {
  return poisson_binomial_lcdf<T_y, T_theta>(y, theta);
}

}  // namespace math
}  // namespace stan
#endif
