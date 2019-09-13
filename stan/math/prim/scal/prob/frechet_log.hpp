#ifndef STAN_MATH_PRIM_SCAL_PROB_FRECHET_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_FRECHET_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/frechet_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>frechet_lpdf</code>
 */
template <bool propto, typename T_y, typename T_shape, typename T_scale>
inline auto frechet_log(T_y&& y, T_shape&& alpha, T_scale&& sigma) {
  return frechet_lpdf<propto>(y, alpha, sigma);
}

/**
 * @deprecated use <code>frechet_lpdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
inline auto frechet_log(T_y&& y, T_shape&& alpha, T_scale&& sigma) {
  return frechet_lpdf(y, alpha, sigma);
}

}  // namespace math
}  // namespace stan
#endif
