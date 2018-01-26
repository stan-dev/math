#ifndef STAN_MATH_PRIM_MAT_PROB_NORMAL_ID_GLM_LOG_HPP
#define STAN_MATH_PRIM_MAT_PROB_NORMAL_ID_GLM_LOG_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/prob/normal_id_glm_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>normal_id_glm_lpdf</code>
 */
template <bool propto, typename T_n, typename T_x, typename T_beta,
          typename T_alpha, typename T_scale>
typename return_type<T_n, T_x, T_beta, T_alpha, T_scale>::type
normal_id_glm_log(const T_n &n, const T_x &x, const T_beta &beta,
                  const T_alpha &alpha, const T_scale &sigma) {
  return normal_id_glm_lpdf<propto, T_n, T_x, T_beta, T_alpha, T_scale>(
      n, x, beta, alpha, sigma);
}

/**
 * @deprecated use <code>normal_id_glm_lpdf</code>
 */
template <typename T_n, typename T_x, typename T_beta, typename T_alpha,
          typename T_scale>
inline typename return_type<T_n, T_x, T_beta, T_alpha, T_scale>::type
normal_id_glm_log(const T_n &n, const T_x &x, const T_beta &beta,
                  const T_alpha &alpha, const T_scale &sigma) {
  return normal_id_glm_lpdf<false>(n, x, beta, alpha, sigma);
}
}  // namespace math
}  // namespace stan
#endif
