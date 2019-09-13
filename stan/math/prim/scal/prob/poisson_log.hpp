#ifndef STAN_MATH_PRIM_SCAL_PROB_POISSON_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_POISSON_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/poisson_lpmf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>poisson_lpmf</code>
 */
template <bool propto, typename T_n, typename T_rate>
inline auto poisson_log(T_n&& n, T_rate&& lambda) {
  return poisson_lpmf<propto>(n, lambda);
}

/**
 * @deprecated use <code>poisson_lpmf</code>
 */
template <typename T_n, typename T_rate>
inline auto poisson_log(T_n&& n, T_rate&& lambda) {
  return poisson_lpmf(n, lambda);
}

}  // namespace math
}  // namespace stan
#endif
