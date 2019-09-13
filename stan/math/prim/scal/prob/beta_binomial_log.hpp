#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_BINOMIAL_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_BINOMIAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/beta_binomial_lpmf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>beta_binomial_lpmf</code>
 */
template <bool propto, typename T_n, typename T_N, typename T_size1,
          typename T_size2>
inline auto beta_binomial_log(T_n&& n, T_N&& N, T_size1&& alpha,
                              T_size2&& beta) {
  return beta_binomial_lpmf<propto>(n, N, alpha, beta);
}

/**
 * @deprecated use <code>beta_binomial_lpmf</code>
 */
template <typename T_n, typename T_N, typename T_size1, typename T_size2>
inline auto beta_binomial_log(T_n&& n, T_N&& N, T_size1&& alpha,
                              T_size2&& beta) {
  return beta_binomial_lpmf(n, N, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
