
#ifndef STAN_MATH_PRIM_PROB_GAMMA_POISSON_LOG_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_POISSON_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/gamma_poisson_lpmf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>gamma_poisson_lpmf</code>
 */
template <bool propto, typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> gamma_poisson_log(const T_n& n,
                                                      const T_shape& alpha,
                                                      const T_inv_scale& beta) {
  return gamma_poisson_lpmf<propto, T_n, T_shape, T_inv_scale>(n, alpha, beta);
}

/** \ingroup prob_dists
 * @deprecated use <code>gamma_poisson_lpmf</code>
 */
template <typename T_n, typename T_shape, typename T_inv_scale>
inline return_type_t<T_shape, T_inv_scale> gamma_poisson_log(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  return gamma_poisson_lpmf<T_n, T_shape, T_inv_scale>(n, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
