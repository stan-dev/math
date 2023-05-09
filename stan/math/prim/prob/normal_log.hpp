#ifndef STAN_MATH_PRIM_PROB_NORMAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/normal_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * the specified mean(s) and deviation(s). y, mu, or sigma can
 * each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @deprecated use <code>normal_lpdf</code>
 *
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * for the normal distribution.
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 * @tparam T_y Underlying type of scalar in sequence.
 * @tparam T_loc Type of location parameter.
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
inline auto normal_log(const T_y& y, const T_loc& mu, const T_scale& sigma) {
  return normal_lpdf<propto>(y, mu, sigma);
}

/** \ingroup prob_dists
 * @deprecated use <code>normal_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto normal_log(const T_y& y, const T_loc& mu, const T_scale& sigma) {
  return normal_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
