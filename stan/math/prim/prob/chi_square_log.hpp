#ifndef STAN_MATH_PRIM_PROB_CHI_SQUARE_LOG_HPP
#define STAN_MATH_PRIM_PROB_CHI_SQUARE_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/chi_square_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>chi_square_lpdf</code>
 */
template <bool propto, typename T_y, typename T_dof>
return_type_t<T_y, T_dof> chi_square_log(const T_y& y, const T_dof& nu) {
  return chi_square_lpdf<propto, T_y, T_dof>(y, nu);
}

/** \ingroup prob_dists
 * @deprecated use <code>chi_square_lpdf</code>
 */
template <typename T_y, typename T_dof>
inline return_type_t<T_y, T_dof> chi_square_log(const T_y& y, const T_dof& nu) {
  return chi_square_lpdf<T_y, T_dof>(y, nu);
}

}  // namespace math
}  // namespace stan
#endif
