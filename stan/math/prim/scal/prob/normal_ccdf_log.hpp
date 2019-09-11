#ifndef STAN_MATH_PRIM_SCAL_PROB_NORMAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NORMAL_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/normal_lccdf.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * @deprecated use <code>normal_lccdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto normal_ccdf_log(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return normal_lccdf(std::forward<T_y>(y), std::forward<T_loc>(mu),
   std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
