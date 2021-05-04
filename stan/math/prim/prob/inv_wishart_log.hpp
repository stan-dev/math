#ifndef STAN_MATH_PRIM_PROB_INV_WISHART_LOG_HPP
#define STAN_MATH_PRIM_PROB_INV_WISHART_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/inv_wishart_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the Inverse-Wishart density for the given W, degrees
 * of freedom, and scale matrix.
 *
 * The scale matrix, S, must be k x k, symmetric, and semi-positive
 * definite.
 *
 * @deprecated use <code>inverse_wishart_lpdf</code>
 *
 * @param W A scalar matrix
 * @param nu Degrees of freedom
 * @param S The scale matrix
 * @return The log of the Inverse-Wishart density at W given nu and S.
 * @throw std::domain_error if nu is not greater than k-1
 * @throw std::domain_error if S is not square, not symmetric, or not
 * semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_dof Type of degrees of freedom.
 * @tparam T_scale Type of scale.
 */
template <bool propto, typename T_y, typename T_dof, typename T_scale>
return_type_t<T_y, T_dof, T_scale> inv_wishart_log(const T_y& W,
                                                   const T_dof& nu,
                                                   const T_scale& S) {
  return inv_wishart_lpdf<propto>(W, nu, S);
}

/** \ingroup multivar_dists
 * @deprecated use <code>inverse_wishart_lpdf</code>
 */
template <typename T_y, typename T_dof, typename T_scale>
inline return_type_t<T_y, T_dof, T_scale> inv_wishart_log(const T_y& W,
                                                          const T_dof& nu,
                                                          const T_scale& S) {
  return inv_wishart_lpdf<>(W, nu, S);
}

}  // namespace math
}  // namespace stan
#endif
