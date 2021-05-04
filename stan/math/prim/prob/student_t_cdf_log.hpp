#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/student_t_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>student_t_lcdf</code>
 */
template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
return_type_t<T_y, T_dof, T_loc, T_scale> student_t_cdf_log(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma) {
  return student_t_lcdf<T_y, T_dof, T_loc, T_scale>(y, nu, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
