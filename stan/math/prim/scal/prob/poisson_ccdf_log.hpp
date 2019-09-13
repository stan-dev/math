#ifndef STAN_MATH_PRIM_SCAL_PROB_POISSON_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_POISSON_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/poisson_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>poisson_lccdf</code>
 */
template <typename T_n, typename T_rate>
inline auto poisson_ccdf_log(T_n&& n, T_rate&& lambda) {
  return poisson_lccdf(n, lambda);
}

}  // namespace math
}  // namespace stan
#endif
