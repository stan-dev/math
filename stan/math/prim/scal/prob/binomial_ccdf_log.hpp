#ifndef STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/binomial_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>binomial_lccdf</code>
 */
template <typename T_n, typename T_N, typename T_prob>
inline auto binomial_ccdf_log(const T_n& n, const T_N& N, const T_prob& theta) {
  return binomial_lccdf(n, N, theta);
}

}  // namespace math
}  // namespace stan
#endif
