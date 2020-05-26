#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/categorical_lpmf.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>categorical_lpmf</code>
 */
template <bool propto, typename T_n, typename T_prob>
inline return_type_t<T_prob> categorical_log(const T_n& ns,
                                             const T_prob& theta) {
  return categorical_lpmf<propto>(ns, theta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>categorical_lpmf</code>
 */
template <typename T_n, typename T_prob>
inline return_type_t<T_prob> categorical_log(const T_n& ns,
                                             const T_prob& theta) {
  return categorical_lpmf(ns, theta);
}

}  // namespace math
}  // namespace stan
#endif
