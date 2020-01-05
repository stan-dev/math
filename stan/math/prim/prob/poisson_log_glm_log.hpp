#ifndef STAN_MATH_PRIM_PROB_POISSON_LOG_GLM_LOG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LOG_GLM_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/poisson_log_glm_lpmf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>poisson_logit_glm_lpmf</code>
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta>
return_type_t<T_x, T_alpha, T_beta> poisson_log_glm_log(const T_y &y,
                                                        const T_x &x,
                                                        const T_alpha &alpha,
                                                        const T_beta &beta) {
  return poisson_log_glm_lpmf<propto, T_y, T_x, T_alpha, T_beta>(y, x, alpha,
                                                                 beta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>poisson_logit_glm_lpmf</code>
 */
template <typename T_y, typename T_x, typename T_alpha, typename T_beta>
inline return_type_t<T_x, T_alpha, T_beta> poisson_log_glm_log(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta) {
  return poisson_log_glm_lpmf<false>(y, x, alpha, beta);
}
}  // namespace math
}  // namespace stan
#endif
