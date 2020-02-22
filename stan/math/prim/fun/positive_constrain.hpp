#ifndef STAN_MATH_PRIM_FUN_POSITIVE_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_POSITIVE_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the positive value for the specified unconstrained input.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x)\f$.
 *
 * @param x Arbitrary input.
 * @return Input transformed to be positive.
 */
template <typename T>
inline auto positive_constrain(T&& x) {
  using std::exp;
  return exp(std::forward<T>(x));
}

/**
 * Return the positive value for the specified unconstrained scalar input,
 * incrementing the scalar reference with the log absolute
 * Jacobian determinant.
 *
 * <p>See <code>positive_constrain(T)</code> for details
 * of the transform.  The log absolute Jacobian determinant is
 *
 * <p>\f$\log | \frac{d}{dx} \mbox{exp}(x) |
 *    = \log | \mbox{exp}(x) | =  x\f$.
 *
 * @tparam T1 type of unconstrained value
 * @tparam T2 type of log prob
 * @param x unconstrained value
 * @param lp log density reference.
 * @return positive constrained version of unconstrained value
 */
template <typename T1, typename T2,
          typename = require_all_stan_scalar_t<T1, T2>>
inline auto positive_constrain(T1&& x, T2&& lp) {
  using std::exp;
  lp += x;
  return exp(std::forward<T1>(x));
}

/**
 * Return the positive value for the specified unconstrained eigen input,
 * incrementing the scalar reference with the log absolute
 * Jacobian determinant.
 *
 * <p>See <code>positive_constrain(T)</code> for details
 * of the transform.  The log absolute Jacobian determinant is
 *
 * <p>\f$\log | \frac{d}{dx} \mbox{exp}(x) |
 *    = \log | \mbox{exp}(x) | =  x\f$.
 *
 * @tparam T1 type of Eigen matrix.
 * @tparam T2 scalar type for log prob.
 * @param x unconstrained value
 * @param lp log density reference.
 * @return positive constrained version of unconstrained value
 */
template <typename T1, typename T2, typename = require_stan_scalar_t<T2>,
          typename = require_eigen_t<T1>>
inline auto positive_constrain(T1&& x, T2&& lp) {
  using std::exp;
  lp += x.sum();
  return exp(std::forward<T1>(x));
}

}  // namespace math
}  // namespace stan

#endif
