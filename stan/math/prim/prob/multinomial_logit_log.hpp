#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_LOG_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/multinomial_logit_lpmf.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>multinomial_logit_lpmf</code>
 */
template <bool propto, typename T_beta, typename T_prob = scalar_type_t<T_beta>,
          require_eigen_col_vector_t<T_beta>* = nullptr>
return_type_t<T_prob> multinomial_logit_log(const std::vector<int>& ns,
                                            const T_beta& beta) {
  return multinomial_logit_lpmf<propto, T_beta>(ns, beta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>multinomial_logit_lpmf</code>
 */
template <typename T_beta, typename T_prob = scalar_type_t<T_beta>,
          require_eigen_col_vector_t<T_beta>* = nullptr>
return_type_t<T_prob> multinomial_logit_log(const std::vector<int>& ns,
                                            const T_beta& beta) {
  return multinomial_logit_lpmf<false>(ns, beta);
}

}  // namespace math
}  // namespace stan
#endif
