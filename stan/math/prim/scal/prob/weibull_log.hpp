#ifndef STAN_MATH_PRIM_SCAL_PROB_WEIBULL_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_WEIBULL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/weibull_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>weibull_lpdf</code>
 */
template <bool propto, typename T_y, typename T_shape, typename T_scale>
inline auto weibull_log(T_y&& y, T_shape&& alpha,
                        T_scale&& sigma) {
  return weibull_lpdf<propto>(std::forward<T_y>(y), std::forward<T_shape>(alpha), std::forward<T_scale>(sigma));
}

/**
 * @deprecated use <code>weibull_lpdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
inline auto weibull_log(T_y&& y, T_shape&& alpha,
                        T_scale&& sigma) {
  return weibull_lpdf(std::forward<T_y>(y), std::forward<T_shape>(alpha), std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
