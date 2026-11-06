#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/lognormal_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>lognormal_lcdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> lognormal_cdf_log(T_y&& y, T_loc&& mu,
                                                            T_scale&& sigma) {
  return lognormal_lcdf(std::forward<T_y>(y), std::forward<T_loc>(mu),
                        std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
