#ifndef STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_LOG_HPP
#define STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_LOG_HPP


#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/prob/bernoulli_logit_glm_lpmf.hpp>

namespace stan {
  namespace math {

    /**
     * @deprecated use <code>bernoulli_logit_glm_lpmf</code>
     */
    template <bool propto, typename T_n, typename T_x, typename T_beta,
              typename T_alpha>
    typename return_type<T_x, T_beta, T_alpha>::type
    bernoulli_logit_glm_log(const T_n &n, const T_x &x, const T_beta &beta,
                             const T_alpha &alpha) {
      return bernoulli_logit_glm_lpmf<propto, T_n, T_x, T_beta, T_alpha>
        (n, x, beta, alpha);
    }

    /**
     * @deprecated use <code>bernoulli_logit_glm_lpmf</code>
     */
    template <typename T_n, typename T_x, typename T_beta, typename T_alpha>
    inline
        typename return_type<T_x, T_beta, T_alpha>::type
        bernoulli_logit_glm_log(const T_n &n, const T_x &x, const T_beta &beta,
                                 const T_alpha &alpha) {
      return bernoulli_logit_glm_lpmf<false>(n, x, beta, alpha);
    }
  }  // namespace math
}  // namespace stan
#endif
