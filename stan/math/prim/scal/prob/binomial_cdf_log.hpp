#ifndef STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/binomial_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>binomial_lcdf</code>
 */
template <typename T_n, typename T_N, typename T_prob>
inline auto binomial_cdf_log(T_n&& n, T_N&& N, T_prob&& theta) {
  return binomial_lcdf(n, N, theta);
}

}  // namespace math
}  // namespace stan
#endif
