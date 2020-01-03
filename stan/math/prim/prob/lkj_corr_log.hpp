#ifndef STAN_MATH_PRIM_PROB_LKJ_CORR_LOG_HPP
#define STAN_MATH_PRIM_PROB_LKJ_CORR_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/prob/lkj_corr_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>lkj_corr_lpdf</code>
 */
template <bool propto, typename T_y, typename T_shape>
return_type_t<T_y, T_shape> lkj_corr_log(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y,
    const T_shape& eta) {
  return lkj_corr_lpdf<propto, T_y, T_shape>(y, eta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>lkj_corr_lpdf</code>
 */
template <typename T_y, typename T_shape>
inline return_type_t<T_y, T_shape> lkj_corr_log(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y,
    const T_shape& eta) {
  return lkj_corr_lpdf<T_y, T_shape>(y, eta);
}

}  // namespace math
}  // namespace stan
#endif
