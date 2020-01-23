#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/std_normal_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of a standard normal density for the specified scalar(s).
 * y can be either a scalar or a vector.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation.
 *
 * @deprecated use <code>std_normal_lpdf</code>
 *
 * @tparam T_y Underlying type of scalar in sequence.
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of the densities.
 * @throw std::domain_error if any scalar is nan.
 */
template <bool propto, typename T_y>
return_type_t<T_y> std_normal_log(const T_y& y) {
  return std_normal_lpdf<propto, T_y>(y);
}

/** \ingroup prob_dists
 * @deprecated use <code>std_normal_lpdf</code>
 */
template <typename T_y>
inline return_type_t<T_y> std_normal_log(const T_y& y) {
  return std_normal_lpdf<T_y>(y);
}

}  // namespace math
}  // namespace stan
#endif
