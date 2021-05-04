#ifndef STAN_MATH_PRIM_PROB_NORMAL_ID_GLM_LOG_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_ID_GLM_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/normal_id_glm_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>normal_id_glm_lpdf</code>
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, typename T_scale>
return_type_t<T_y, T_x, T_alpha, T_beta, T_scale> normal_id_glm_log(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta,
    const T_scale &sigma) {
  return normal_id_glm_lpdf<propto, T_y, T_x, T_alpha, T_beta, T_scale>(
      y, x, alpha, beta, sigma);
}

/** \ingroup multivar_dists
 * @deprecated use <code>normal_id_glm_lpdf</code>
 */
template <typename T_y, typename T_x, typename T_alpha, typename T_beta,
          typename T_scale>
inline return_type_t<T_y, T_x, T_alpha, T_beta, T_scale> normal_id_glm_log(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta,
    const T_scale &sigma) {
  return normal_id_glm_lpdf<false>(y, x, alpha, beta, sigma);
}
}  // namespace math
}  // namespace stan
#endif
