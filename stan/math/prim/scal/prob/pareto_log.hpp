#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/pareto_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>pareto_lpdf</code>
 */
template <bool propto, typename T_y, typename T_scale, typename T_shape>
inline auto pareto_log(T_y&& y, T_scale&& y_min,
                       T_shape&& alpha) {
  return pareto_lpdf<propto>(y, y_min, alpha);
}

/**
 * @deprecated use <code>pareto_lpdf</code>
 */
template <typename T_y, typename T_scale, typename T_shape>
inline auto pareto_log(T_y&& y, T_scale&& y_min,
                       T_shape&& alpha) {
  return pareto_lpdf(y, y_min, alpha);
}

}  // namespace math
}  // namespace stan
#endif
