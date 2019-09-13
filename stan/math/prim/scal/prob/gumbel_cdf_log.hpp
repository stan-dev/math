#ifndef STAN_MATH_PRIM_SCAL_PROB_GUMBEL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_GUMBEL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/gumbel_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>gumbel_lcdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto gumbel_cdf_log(T_y&& y, T_loc&& mu, T_scale&& beta) {
  return gumbel_lcdf(y, mu, beta);
}

}  // namespace math
}  // namespace stan
#endif
