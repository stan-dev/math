#ifndef STAN_MATH_PRIM_SCAL_PROB_STUDENT_T_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_STUDENT_T_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/student_t_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>student_t_lccdf</code>
 */
template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
inline auto student_t_ccdf_log(const T_y& y, const T_dof& nu, const T_loc& mu,
                               const T_scale& sigma) {
  return student_t_lccdf(y, nu, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
