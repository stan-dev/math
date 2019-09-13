#ifndef STAN_MATH_PRIM_SCAL_PROB_STUDENT_T_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_STUDENT_T_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/student_t_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>student_t_lcdf</code>
 */
template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
inline auto student_t_cdf_log(T_y&& y, T_dof&& nu, T_loc&& mu,
                              T_scale&& sigma) {
  return student_t_lcdf(y, nu, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
