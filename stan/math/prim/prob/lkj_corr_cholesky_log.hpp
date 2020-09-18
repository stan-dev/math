#ifndef STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_LOG_HPP
#define STAN_MATH_PRIM_PROB_LKJ_CORR_CHOLESKY_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/lkj_corr_cholesky_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>lkj_corr_cholesky_lpdf</code>
 */
template <bool propto, typename T_covar, typename T_shape>
return_type_t<T_covar, T_shape> lkj_corr_cholesky_log(const T_covar& L,
                                                      const T_shape& eta) {
  return lkj_corr_cholesky_lpdf<propto>(L, eta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>lkj_corr_cholesky_lpdf</code>
 */
template <typename T_covar, typename T_shape>
inline return_type_t<T_covar, T_shape> lkj_corr_cholesky_log(
    const T_covar& L, const T_shape& eta) {
  return lkj_corr_cholesky_lpdf<>(L, eta);
}

}  // namespace math
}  // namespace stan
#endif
