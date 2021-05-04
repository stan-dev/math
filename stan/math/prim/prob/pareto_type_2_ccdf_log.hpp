#ifndef STAN_MATH_PRIM_PROB_PARETO_TYPE_2_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_PARETO_TYPE_2_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/pareto_type_2_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>pareto_type_2_lccdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_ccdf_log(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  return pareto_type_2_lccdf<T_y, T_loc, T_scale, T_shape>(y, mu, lambda,
                                                           alpha);
}

}  // namespace math
}  // namespace stan
#endif
