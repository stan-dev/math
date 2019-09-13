#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/neg_binomial_lpmf.hpp>
#include <utility>
namespace stan {
namespace math {

/**
 * @deprecated use <code>neg_binomial_lpmf</code>
 */
template <bool propto, typename T_n, typename T_shape, typename T_inv_scale>
inline auto neg_binomial_log(T_n&& n, T_shape&& alpha, T_inv_scale&& beta) {
  return neg_binomial_lpmf<propto>(std::forward<T_n>(n),
                                   std::forward<T_shape>(alpha),
                                   std::forward<T_inv_scale>(beta));
}

/**
 * @deprecated use <code>neg_binomial_lpmf</code>
 */
template <typename T_n, typename T_shape, typename T_inv_scale>
inline auto neg_binomial_log(T_n&& n, T_shape&& alpha, T_inv_scale&& beta) {
  return neg_binomial_lpmf(std::forward<T_n>(n), std::forward<T_shape>(alpha),
                           std::forward<T_inv_scale>(beta));
}

}  // namespace math
}  // namespace stan
#endif
