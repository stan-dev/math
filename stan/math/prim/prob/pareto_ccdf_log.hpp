#ifndef STAN_MATH_PRIM_PROB_PARETO_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_PARETO_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/pareto_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>pareto_lccdf</code>
 */
template <typename T_y, typename T_scale, typename T_shape>
return_type_t<T_y, T_scale, T_shape> pareto_ccdf_log(const T_y& y,
                                                     const T_scale& y_min,
                                                     const T_shape& alpha) {
  return pareto_lccdf<T_y, T_scale, T_shape>(y, y_min, alpha);
}

}  // namespace math
}  // namespace stan
#endif
