#ifndef STAN_MATH_PRIM_PROB_POISSON_GAMMA_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_GAMMA_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/poisson_gamma_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>poisson_gamma_lccdf</code>
 */
template <typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> poisson_gamma_ccdf_log(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  return poisson_gamma_lccdf<T_n, T_shape, T_inv_scale>(n, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
