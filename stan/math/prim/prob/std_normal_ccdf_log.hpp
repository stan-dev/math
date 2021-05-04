#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/std_normal_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>std_normal_lccdf</code>
 */
template <typename T_y>
inline return_type_t<T_y> std_normal_ccdf_log(const T_y& y) {
  return std_normal_lccdf<T_y>(y);
}

}  // namespace math
}  // namespace stan
#endif
