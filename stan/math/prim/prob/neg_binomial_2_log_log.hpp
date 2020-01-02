#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_LOG_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/neg_binomial_2_log_lpmf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>neg_binomial_2_log_lpmf</code>
 */
template <bool propto, typename T_n, typename T_log_location,
          typename T_precision>
return_type_t<T_log_location, T_precision> neg_binomial_2_log_log(
    const T_n& n, const T_log_location& eta, const T_precision& phi) {
  return neg_binomial_2_log_lpmf<propto, T_n, T_log_location, T_precision>(
      n, eta, phi);
}

/** \ingroup prob_dists
 * @deprecated use <code>neg_binomial_2_log_lpmf</code>
 */
template <typename T_n, typename T_log_location, typename T_precision>
inline return_type_t<T_log_location, T_precision> neg_binomial_2_log_log(
    const T_n& n, const T_log_location& eta, const T_precision& phi) {
  return neg_binomial_2_log_lpmf<T_n, T_log_location, T_precision>(n, eta, phi);
}
}  // namespace math
}  // namespace stan
#endif
