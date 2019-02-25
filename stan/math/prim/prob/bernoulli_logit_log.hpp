#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_LOG_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_LOG_HPP

#include <stanh/prim/meta/return_type.hpp>
#include <stanh/prim/prob/bernoulli_logit_lpmf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>bernoulli_logit_lpmf</code>
 */
template <bool propto, typename T_n, typename T_prob>
typename return_type<T_prob>::type bernoulli_logit_log(const T_n& n,
                                                       const T_prob& theta) {
  return bernoulli_logit_lpmf<propto, T_n, T_prob>(n, theta);
}

/**
 * @deprecated use <code>bernoulli_logit_lpmf</code>
 */
template <typename T_n, typename T_prob>
inline typename return_type<T_prob>::type bernoulli_logit_log(
    const T_n& n, const T_prob& theta) {
  return bernoulli_logit_lpmf<T_n, T_prob>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
