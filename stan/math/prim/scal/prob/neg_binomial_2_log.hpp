#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/neg_binomial_2_lpmf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>neg_binomial_2_lpmf</code>
 */
template <bool propto, typename T_n, typename T_location, typename T_precision>
inline auto neg_binomial_2_log(T_n&& n, T_location&& mu,
                               T_precision&& phi) {
  return neg_binomial_2_lpmf<propto>(n, mu, phi);
}

/**
 * @deprecated use <code>neg_binomial_2_lpmf</code>
 */
template <typename T_n, typename T_location, typename T_precision>
inline auto neg_binomial_2_log(T_n&& n, T_location&& mu,
                               T_precision&& phi) {
  return neg_binomial_2_lpmf(n, mu, phi);
}

}  // namespace math
}  // namespace stan
#endif
