#ifndef STAN_MATH_PRIM_SCAL_PROB_RAYLEIGH_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_RAYLEIGH_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/rayleigh_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>rayleigh_lccdf</code>
 */
template <typename T_y, typename T_scale>
inline auto rayleigh_ccdf_log(T_y&& y, T_scale&& sigma) {
  return rayleigh_lccdf(y, sigma);
}

}  // namespace math
}  // namespace stan
#endif
