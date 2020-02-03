#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/neg_binomial_2_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>neg_binomial_2_lccdf</code>
 */
template <typename T_n, typename T_location, typename T_precision>
return_type_t<T_location, T_precision> neg_binomial_2_ccdf_log(
    const T_n& n, const T_location& mu, const T_precision& phi) {
  return neg_binomial_2_lccdf<T_n, T_location, T_precision>(n, mu, phi);
}

}  // namespace math
}  // namespace stan
#endif
