#ifndef STAN_MATH_PRIM_SCAL_PROB_WEIBULL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_WEIBULL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/weibull_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>weibull_lcdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
inline auto weibull_cdf_log(T_y&& y, T_shape&& alpha, T_scale&& sigma) {
  return weibull_lcdf(std::forward<T_y>(y), std::forward<T_shape>(alpha),
                      std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
