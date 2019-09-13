#ifndef STAN_MATH_PRIM_SCAL_PROB_NORMAL_SUFFICIENT_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NORMAL_SUFFICIENT_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/normal_sufficient_lpdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>normal_lpdf</code>
 */
template <bool propto, typename T_y, typename T_s, typename T_n, typename T_loc,
          typename T_scale>
inline auto normal_sufficient_log(T_y&& y_bar, T_s&& s_squared,
                                  T_n&& n_obs, T_loc&& mu,
                                  T_scale&& sigma) {
  return normal_sufficient_lpdf<propto>(y_bar, s_squared, n_obs, mu, sigma);
}

/**
 * @deprecated use <code>normal_lpdf</code>
 */
template <typename T_y, typename T_s, typename T_n, typename T_loc,
          typename T_scale>
inline auto normal_sufficient_log(T_y&& y_bar, T_s&& s_squared,
                                  T_n&& n_obs, T_loc&& mu,
                                  T_scale&& sigma) {
  return normal_sufficient_lpdf(y_bar, s_squared, n_obs, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
