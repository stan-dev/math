#ifndef STAN_MATH_PRIM_SCAL_PROB_SCALED_INV_CHI_SQUARE_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_SCALED_INV_CHI_SQUARE_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/scaled_inv_chi_square_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>scaled_inv_chi_square_lcdf</code>
 */
template <typename T_y, typename T_dof, typename T_scale>
inline auto scaled_inv_chi_square_cdf_log(T_y&& y, T_dof&& nu,
                                          T_scale&& s) {
  return scaled_inv_chi_square_lcdf(y, nu, s);
}

}  // namespace math
}  // namespace stan
#endif
