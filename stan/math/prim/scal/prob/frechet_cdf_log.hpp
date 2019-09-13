#ifndef STAN_MATH_PRIM_SCAL_PROB_FRECHET_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_FRECHET_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/frechet_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>frechet_lcdf</code>
 */
template <typename T_y, typename T_shape, typename T_scale>
inline auto frechet_cdf_log(T_y&& y, T_shape&& alpha,
                            T_scale&& sigma) {
  return frechet_lcdf(y, alpha, sigma);
}

}  // namespace math
}  // namespace stan
#endif
