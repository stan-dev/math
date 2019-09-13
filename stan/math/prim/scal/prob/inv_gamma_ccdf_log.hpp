#ifndef STAN_MATH_PRIM_SCAL_PROB_INV_GAMMA_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_INV_GAMMA_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/inv_gamma_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>inv_gamma_lccdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
inline auto inv_gamma_ccdf_log(T_y&& y, T_shape&& alpha,
                               T_scale&& beta) {
  return inv_gamma_lccdf(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
