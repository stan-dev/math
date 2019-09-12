#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/bernoulli_lpmf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>bernoulli_lpmf</code>
 */
template <bool propto, typename T_n, typename T_prob>
inline auto bernoulli_log(const T_n& n, const T_prob& theta) {
  return bernoulli_lpmf<propto>(n, theta);
}

/**
 * @deprecated use <code>bernoulli_lpmf</code>
 */
template <typename T_y, typename T_prob>
inline auto bernoulli_log(const T_y& n, const T_prob& theta) {
  return bernoulli_lpmf(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
