#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/bernoulli_lpmf.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * @deprecated use <code>bernoulli_lpmf</code>
 */
template <bool propto, typename T_n, typename T_prob>
inline auto bernoulli_log(T_n&& n, T_prob&& theta) {
  return bernoulli_lpmf<propto>(std::forward<T_n>(n),
                                std::forward<T_prob>(theta));
}

/**
 * @deprecated use <code>bernoulli_lpmf</code>
 */
template <typename T_n, typename T_prob>
inline auto bernoulli_log(T_n&& n, T_prob&& theta) {
  return bernoulli_lpmf(std::forward<T_n>(n), std::forward<T_prob>(theta));
}

}  // namespace math
}  // namespace stan
#endif
