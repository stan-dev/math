#ifndef STAN_MATH_PRIM_FUN_POSITIVE_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_POSITIVE_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
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
 * @param x Arbitrary input scalar or container.
 * @return Input transformed to be positive.
 */
template <typename T>
inline auto positive_constrain(const T& x) {
  return exp(x);
}

/**
 * Return the positive value for the specified unconstrained input,
 * incrementing the scalar reference with the log absolute
 * Jacobian determinant.
 *
 * <p>See <code>positive_constrain(T)</code> for details
 * of the transform.  The log absolute Jacobian determinant is
 *
 * <p>\f$\log | \frac{d}{dx} \mbox{exp}(x) |
 *    = \log | \mbox{exp}(x) | =  x\f$.
 *
 * @tparam T type of unconstrained value
 * @param x unconstrained value or container
 * @param lp log density reference.
 * @return positive constrained version of unconstrained value(s)
 */
template <typename T, typename S>
inline auto positive_constrain(const T& x, S& lp) {
  lp += sum(x);
  return exp(x);
}

}  // namespace math
}  // namespace stan

#endif
