#ifndef STAN_MATH_PRIM_PROB_GAMMA_POISSON_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_POISSON_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/gamma_poisson_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>gamma_poisson_lcdf</code>
 */
template <typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> gamma_poisson_cdf_log(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  return gamma_poisson_lcdf<T_n, T_shape, T_inv_scale>(n, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
