#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/neg_binomial_2_log_lpmf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>neg_binomial_2_log_lpmf</code>
 */
template <bool propto, typename T_n, typename T_log_location,
          typename T_precision>
inline auto neg_binomial_2_log_log(T_n&& n, T_log_location&& eta,
                                   T_precision&& phi) {
  return neg_binomial_2_log_lpmf<propto>(n, eta, phi);
}

/**
 * @deprecated use <code>neg_binomial_2_log_lpmf</code>
 */
template <typename T_n, typename T_log_location, typename T_precision>
inline auto neg_binomial_2_log_log(T_n&& n, T_log_location&& eta,
                                   T_precision&& phi) {
  return neg_binomial_2_log_lpmf(n, eta, phi);
}
}  // namespace math
}  // namespace stan
#endif
