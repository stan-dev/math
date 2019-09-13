#ifndef STAN_MATH_PRIM_SCAL_PROB_RAYLEIGH_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_RAYLEIGH_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/rayleigh_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>rayleigh_lpdf</code>
 */
template <bool propto, typename T_y, typename T_scale>
inline auto rayleigh_log(const T_y& y, const T_scale& sigma) {
  return rayleigh_lpdf<propto>(y, sigma);
}

/**
 * @deprecated use <code>rayleigh_lpdf</code>
 */
template <typename T_y, typename T_scale>
inline auto rayleigh_log(const T_y& y, const T_scale& sigma) {
  return rayleigh_lpdf(y, sigma);
}

}  // namespace math
}  // namespace stan
#endif
