#ifndef STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/lognormal_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>lognormal_lpdf</code>
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
inline auto lognormal_log(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return lognormal_lpdf<propto>(std::forward<T_y>(y), std::forward<T_loc>(mu),
                                std::forward<T_scale>(sigma));
}

/**
 * @deprecated use <code>lognormal_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto lognormal_log(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return lognormal_lpdf(std::forward<T_y>(y), std::forward<T_loc>(mu),
                        std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
