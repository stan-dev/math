#ifndef STAN_MATH_PRIM_SCAL_PROB_CAUCHY_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CAUCHY_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/cauchy_lccdf.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * @deprecated use <code>cauchy_lccdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto cauchy_ccdf_log(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  return cauchy_lccdf(std::forward<T_y>(y), std::forward<T_loc>(mu),
                      std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
