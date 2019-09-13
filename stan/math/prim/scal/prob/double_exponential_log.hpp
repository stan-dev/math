#ifndef STAN_MATH_PRIM_SCAL_PROB_DOUBLE_EXPONENTIAL_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_DOUBLE_EXPONENTIAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/double_exponential_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>double_exponential_lpdf</code>
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
inline auto double_exponential_log(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return double_exponential_lpdf<propto>(y, mu, sigma);
}

/**
 * @deprecated use <code>double_exponential_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto double_exponential_log(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return double_exponential_lpdf(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
