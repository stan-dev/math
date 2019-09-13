#ifndef STAN_MATH_PRIM_SCAL_PROB_GAMMA_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_GAMMA_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/gamma_lccdf.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * @deprecated use <code>gamma_lccdf</code>
 */
template <typename T_y, typename T_shape, typename T_inv_scale>
inline auto gamma_ccdf_log(T_y&& y, T_shape&& alpha,
                           T_inv_scale&& beta) {
  return gamma_lccdf(std::forward<T_y>(y), std::forward<T_shape>(alpha), std::forward<T_inv_scale>(beta));
}

}  // namespace math
}  // namespace stan

#endif
