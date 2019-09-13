#ifndef STAN_MATH_PRIM_SCAL_PROB_LOGISTIC_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_LOGISTIC_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/logistic_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>logistic_lcdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto logistic_cdf_log(T_y&& y, T_loc&& mu,
                             T_scale&& sigma) {
  return logistic_lcdf(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
