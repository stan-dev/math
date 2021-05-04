#ifndef STAN_MATH_PRIM_PROB_RAYLEIGH_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_RAYLEIGH_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/rayleigh_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>rayleigh_lccdf</code>
 */
template <typename T_y, typename T_scale>
return_type_t<T_y, T_scale> rayleigh_ccdf_log(const T_y& y,
                                              const T_scale& sigma) {
  return rayleigh_lccdf<T_y, T_scale>(y, sigma);
}

}  // namespace math
}  // namespace stan
#endif
