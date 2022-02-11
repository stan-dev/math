#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_LOG_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/categorical_logit_lpmf.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>categorical_logit_lpmf</code>
 */
template <bool propto, typename T_n, typename T_prob>
inline return_type_t<T_prob> categorical_logit_log(const T_n& ns,
                                                   const T_prob& beta) {
  return categorical_logit_lpmf<propto>(ns, beta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>categorical_logit_lpmf</code>
 */
template <typename T_n, typename T_prob>
inline return_type_t<T_prob> categorical_logit_log(const T_n& ns,
                                                   const T_prob& beta) {
  return categorical_logit_lpmf(ns, beta);
}

}  // namespace math
}  // namespace stan
#endif
