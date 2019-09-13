#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/neg_binomial_lcdf.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>neg_binomial_lcdf</code>
 */
template <typename T_n, typename T_shape, typename T_inv_scale>
inline auto neg_binomial_cdf_log(T_n&& n, T_shape&& alpha, T_inv_scale&& beta) {
  return neg_binomial_lcdf(std::forward<T_n>(n), std::forward<T_shape>(alpha),
                           std::forward<T_inv_scale>(beta));
}

}  // namespace math
}  // namespace stan
#endif
