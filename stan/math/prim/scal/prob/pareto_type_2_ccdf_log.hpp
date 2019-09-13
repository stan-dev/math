#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/pareto_type_2_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>pareto_type_2_lccdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline auto pareto_type_2_ccdf_log(T_y&& y, T_loc&& mu, T_scale&& lambda,
                                   T_shape&& alpha) {
  return pareto_type_2_lccdf(y, mu, lambda, alpha);
}

}  // namespace math
}  // namespace stan
#endif
