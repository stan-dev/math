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
template <bool propto, typename T_prob>
return_type_t<T_prob> categorical_log(
    int n, const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& theta) {
  return categorical_lpmf<propto, T_prob>(n, theta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>categorical_lpmf</code>
 */
template <typename T_prob>
return_type_t<T_prob> categorical_log(
    const typename math::index_type<Eigen::VectorXd>::type n,
    const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& theta) {
  return categorical_lpmf<T_prob>(n, theta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>categorical_lpmf</code>
 */
template <bool propto, typename T_prob>
return_type_t<T_prob> categorical_log(
    const std::vector<int>& ns,
    const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& theta) {
  return categorical_lpmf<propto, T_prob>(ns, theta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>categorical_lpmf</code>
 */
template <typename T_prob>
inline return_type_t<T_prob> categorical_log(
    const std::vector<int>& ns,
    const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& theta) {
  return categorical_lpmf<false>(ns, theta);
}

}  // namespace math
}  // namespace stan
#endif
