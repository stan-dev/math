#ifndef STAN_MATH_PRIM_PROB_GAUSSIAN_DLM_OBS_LOG_HPP
#define STAN_MATH_PRIM_PROB_GAUSSIAN_DLM_OBS_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/gaussian_dlm_obs_lpdf.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {
/** \ingroup multivar_dists
 * @deprecated use <code>gaussian_dlm_obs_lpdf</code>
 */
template <bool propto, typename T_y, typename T_F, typename T_G, typename T_V,
          typename T_W, typename T_m0, typename T_C0>
inline return_type_t<T_y, T_F, T_G, T_V, T_W, T_m0, T_C0> gaussian_dlm_obs_log(
    const T_y& y, const T_F& F, const T_G& G, const T_V& V, const T_W& W,
    const T_m0& m0, const T_C0& C0) {
  return gaussian_dlm_obs_lpdf<propto>(y, F, G, V, W, m0, C0);
}

/** \ingroup multivar_dists
 * @deprecated use <code>gaussian_dlm_obs_lpdf</code>
 */
template <typename T_y, typename T_F, typename T_G, typename T_V, typename T_W,
          typename T_m0, typename T_C0>
inline return_type_t<T_y, T_F, T_G, T_V, T_W, T_m0, T_C0> gaussian_dlm_obs_log(
    const T_y& y, const T_F& F, const T_G& G, const T_V& V, const T_W& W,
    const T_m0& m0, const T_C0& C0) {
  return gaussian_dlm_obs_lpdf<>(y, F, G, V, W, m0, C0);
}

}  // namespace math
}  // namespace stan
#endif
