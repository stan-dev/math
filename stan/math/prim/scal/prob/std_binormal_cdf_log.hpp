#ifndef STAN_MATH_PRIM_SCAL_PROB_STD_BINORMAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_STD_BINORMAL_CDF_LOG_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/prob/std_binormal_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>normal_lcdf</code>
 */
template <typename T_y_1, typename T_y_2, typename T_rho>
typename return_type<T_y_1, T_y_2, T_rho>::type std_binormal_cdf_log(
    const T_y_1& y_1, const T_y_2& y_2, const T_rho& rho) {
  return std_binormal_lcdf<T_y_1, T_y_2, T_rho>(y_1, y_2, rho);
}

}  // namespace math
}  // namespace stan
#endif
