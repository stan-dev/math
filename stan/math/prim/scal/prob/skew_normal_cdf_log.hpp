#ifndef STAN_MATH_PRIM_SCAL_PROB_SKEW_NORMAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_SKEW_NORMAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/skew_normal_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>skew_normal_lcdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline auto skew_normal_cdf_log(T_y&& y, T_loc&& mu, T_scale&& sigma,
                                T_shape&& alpha) {
  return skew_normal_lcdf(y, mu, sigma, alpha);
}

}  // namespace math
}  // namespace stan
#endif
