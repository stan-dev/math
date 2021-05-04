#ifndef STAN_MATH_PRIM_PROB_UNIFORM_LOG_HPP
#define STAN_MATH_PRIM_PROB_UNIFORM_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/uniform_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of a uniform density for the given
 * y, lower, and upper bound.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{U}}(\alpha, \beta) \\
 \log (p (y \, |\, \alpha, \beta)) &=& \log \left( \frac{1}{\beta-\alpha}
 \right) \\
 &=& \log (1) - \log (\beta - \alpha) \\
 &=& -\log (\beta - \alpha) \\
 & & \mathrm{ where } \; y \in [\alpha, \beta], \log(0) \; \mathrm{otherwise}
 \f}
 *
 * @deprecated use <code>uniform_lpdf</code>
 *
 * @param y A scalar variable.
 * @param alpha Lower bound.
 * @param beta Upper bound.
 * @throw std::invalid_argument if the lower bound is greater than
 *    or equal to the lower bound
 * @tparam T_y Type of scalar.
 * @tparam T_low Type of lower bound.
 * @tparam T_high Type of upper bound.
 */
template <bool propto, typename T_y, typename T_low, typename T_high>
return_type_t<T_y, T_low, T_high> uniform_log(const T_y& y, const T_low& alpha,
                                              const T_high& beta) {
  return uniform_lpdf<propto, T_y, T_low, T_high>(y, alpha, beta);
}

/** \ingroup prob_dists
 * @deprecated use <code>uniform_lpdf</code>
 */
template <typename T_y, typename T_low, typename T_high>
inline return_type_t<T_y, T_low, T_high> uniform_log(const T_y& y,
                                                     const T_low& alpha,
                                                     const T_high& beta) {
  return uniform_lpdf<T_y, T_low, T_high>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
