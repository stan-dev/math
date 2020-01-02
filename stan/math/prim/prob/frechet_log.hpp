#ifndef STAN_MATH_PRIM_PROB_FRECHET_LOG_HPP
#define STAN_MATH_PRIM_PROB_FRECHET_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/frechet_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>frechet_lpdf</code>
 */
template <bool propto, typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> frechet_log(const T_y& y,
                                                 const T_shape& alpha,
                                                 const T_scale& sigma) {
  return frechet_lpdf<propto, T_y, T_shape, T_scale>(y, alpha, sigma);
}

/** \ingroup prob_dists
 * @deprecated use <code>frechet_lpdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
inline return_type_t<T_y, T_shape, T_scale> frechet_log(const T_y& y,
                                                        const T_shape& alpha,
                                                        const T_scale& sigma) {
  return frechet_lpdf<T_y, T_shape, T_scale>(y, alpha, sigma);
}

}  // namespace math
}  // namespace stan
#endif
