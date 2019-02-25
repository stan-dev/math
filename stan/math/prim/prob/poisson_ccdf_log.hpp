#ifndef STAN_MATH_PRIM_PROB_POISSON_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_CCDF_LOG_HPP

#include <stanh/prim/meta/return_type.hpp>
#include <stanh/prim/prob/poisson_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>poisson_lccdf</code>
 */
template <typename T_n, typename T_rate>
typename return_type<T_rate>::type poisson_ccdf_log(const T_n& n,
                                                    const T_rate& lambda) {
  return poisson_lccdf<T_n, T_rate>(n, lambda);
}

}  // namespace math
}  // namespace stan
#endif
