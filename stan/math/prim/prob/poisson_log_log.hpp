#ifndef STAN_MATH_PRIM_PROB_POISSON_LOG_LOG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LOG_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/poisson_log_lpmf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>poisson_log_lpmf</code>
 */
template <bool propto, typename T_n, typename T_log_rate>
return_type_t<T_log_rate> poisson_log_log(const T_n& n,
                                          const T_log_rate& alpha) {
  return poisson_log_lpmf<propto, T_n, T_log_rate>(n, alpha);
}

/** \ingroup prob_dists
 * @deprecated use <code>poisson_log_lpmf</code>
 */
template <typename T_n, typename T_log_rate>
inline return_type_t<T_log_rate> poisson_log_log(const T_n& n,
                                                 const T_log_rate& alpha) {
  return poisson_log_lpmf<T_n, T_log_rate>(n, alpha);
}

}  // namespace math
}  // namespace stan
#endif
