#ifndef STAN_MATH_PRIM_SCAL_PROB_NORMAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NORMAL_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/normal_lccdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>normal_lccdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_ccdf_log(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  return normal_lccdf<T_y, T_loc, T_scale>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
