#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_LOG_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/loglogistic_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>loglogistic_lpdf</code>
 */
template <bool propto, typename T_y, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_scale, T_shape> loglogistic_log(const T_y& y,
                                                     const T_scale& alpha,
                                                     const T_shape& beta) {
  return loglogistic_lpdf<propto, T_y, T_scale, T_shape>(y, alpha, beta);
}

/** \ingroup prob_dists
 * @deprecated use <code>loglogistic_lpdf</code>
 */
template <typename T_y, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_scale, T_shape> loglogistic_log(
    const T_y& y, const T_scale& alpha, const T_shape& beta) {
  return loglogistic_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
