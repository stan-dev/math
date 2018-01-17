#ifndef STAN_MATH_PRIM_MAT_PROB_NEG_BINOMIAL_2_LOG_GLM_LOG_HPP
#define STAN_MATH_PRIM_MAT_PROB_NEG_BINOMIAL_2_LOG_GLM_LOG_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/prob/neg_binomial_2_log_glm_lpmf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>neg_binomial_2_log_glm_lpmf</code>
 */
template <bool propto, typename T_n, typename T_x, typename T_beta,
          typename T_alpha, typename T_precision>
typename return_type<T_x, T_beta, T_alpha, T_precision>::type
neg_binomial_2_log_glm_log(const T_n &n, const T_x &x, const T_beta &beta,
                           const T_alpha &alpha, const T_precision &phi) {
  return neg_binomial_2_log_glm_lpmf<propto, T_n, T_x, T_beta, T_alpha,
                                     T_precision>(n, x, beta, alpha, phi);
}

/**
 * @deprecated use <code>poisson_logit_glm_lpmf</code>
 */
template <typename T_n, typename T_x, typename T_beta, typename T_alpha,
          typename T_precision>
inline typename return_type<T_x, T_beta, T_alpha, T_precision>::type
neg_binomial_2_log_glm_log(const T_n &n, const T_x &x, const T_beta &beta,
                           const T_alpha &alpha, const T_precision &phi) {
  return neg_binomial_2_log_glm_lpmf<false>(n, x, beta, alpha, phi);
}
}  // namespace math
}  // namespace stan
#endif
