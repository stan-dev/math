#ifndef STAN_MATH_PRIM_PROB_INV_CHI_SQUARE_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_INV_CHI_SQUARE_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/inv_chi_square_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>inv_chi_square_lccdf</code>
 */
template <typename T_y, typename T_dof>
return_type_t<T_y, T_dof> inv_chi_square_ccdf_log(const T_y& y,
                                                  const T_dof& nu) {
  return inv_chi_square_lccdf<T_y, T_dof>(y, nu);
}

}  // namespace math
}  // namespace stan
#endif
