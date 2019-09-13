#ifndef STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/exponential_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>exponential_lcdf</code>
 */
template <typename T_y, typename T_inv_scale>
inline auto exponential_cdf_log(T_y&& y, T_inv_scale&& beta) {
  return exponential_lcdf(y, beta);
}

}  // namespace math
}  // namespace stan
#endif
