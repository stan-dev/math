#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LOG_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/bernoulli_lpmf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>bernoulli_lpmf</code>
 */
template <bool propto, typename T_n, typename T_prob>
return_type_t<T_prob> bernoulli_log(const T_n& n, const T_prob& theta) {
  return bernoulli_lpmf<propto, T_n, T_prob>(n, theta);
}

/** \ingroup prob_dists
 * @deprecated use <code>bernoulli_lpmf</code>
 */
template <typename T_y, typename T_prob>
inline return_type_t<T_prob> bernoulli_log(const T_y& n, const T_prob& theta) {
  return bernoulli_lpmf<T_y, T_prob>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
