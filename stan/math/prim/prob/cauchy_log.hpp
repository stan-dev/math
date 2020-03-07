#ifndef STAN_MATH_PRIM_PROB_CAUCHY_LOG_HPP
#define STAN_MATH_PRIM_PROB_CAUCHY_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/cauchy_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the Cauchy density for the specified scalar(s) given
 * the specified location parameter(s) and scale parameter(s). y,
 * mu, or sigma can each either be scalar a vector.  Any vector
 * inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/mu/sigma triple.
 *
 * @deprecated use <code>cauchy_lpdf</code>
 *
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 * @tparam T_scale Type of scale.
 * @tparam T_y Type of scalar outcome.
 * @tparam T_loc Type of location.
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> cauchy_log(const T_y& y, const T_loc& mu,
                                              const T_scale& sigma) {
  return cauchy_lpdf<propto, T_y, T_loc, T_scale>(y, mu, sigma);
}

/** \ingroup prob_dists
 * @deprecated use <code>cauchy_lpdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> cauchy_log(const T_y& y,
                                                     const T_loc& mu,
                                                     const T_scale& sigma) {
  return cauchy_lpdf<T_y, T_loc, T_scale>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
