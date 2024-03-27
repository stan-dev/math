#ifndef STAN_MATH_PRIM_CONSTRAINT_PROB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_PROB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log_inv_logit.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log1m_inv_logit.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a probability value constrained to fall between 0 and 1
 * (inclusive) for the specified free scalar.
 *
 * <p>The transform is the inverse logit,
 *
 * <p>\f$f(x) = \mbox{logit}^{-1}(x) = \frac{1}{1 + \exp(x)}\f$.
 *
 * @tparam T type of scalar
 * @param[in] x unconstrained value
 * @return result constrained to fall in (0, 1)
 */
template <typename T>
inline T prob_constrain(const T& x) {
  return inv_logit(x);
}

/**
 * Return a probability value constrained to fall between 0 and 1
 * (inclusive) for the specified free scalar and increment the
 * specified log probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as defined for <code>prob_constrain(T)</code>.
 * The log absolute Jacobian determinant is
 *
 * <p>The log absolute Jacobian determinant is
 *
 * <p>\f$\log | \frac{d}{dx} \mbox{logit}^{-1}(x) |\f$
 * <p>\f$\log ((\mbox{logit}^{-1}(x)) (1 - \mbox{logit}^{-1}(x))\f$
 * <p>\f$\log (\mbox{logit}^{-1}(x)) + \log (1 - \mbox{logit}^{-1}(x))\f$.
 *
 * @tparam T type of scalar
 * @param[in] x unconstrained value
 * @param[in, out] lp log density
 * @return result constrained to fall in (0, 1)
 */
template <typename T>
inline T prob_constrain(const T& x, T& lp) {
  T log_inv_logit_x = log_inv_logit(x);
  lp += log_inv_logit_x + log1m_inv_logit(x);
  return exp(log_inv_logit_x);
}

/**
 * Return a probability value constrained to fall between 0 and 1 (inclusive)
 * for the specified free scalar. If the `Jacobian` parameter is `true`, the log
 * density accumulator is incremented with the log absolute Jacobian determinant
 * of the transform.  All of the transforms are specified with their Jacobians
 * in the *Stan Reference Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T type of scalar
 * @param[in] x unconstrained value
 * @param[in, out] lp log density accumulator
 * @return result constrained to fall in (0, 1)
 */
template <bool Jacobian, typename T>
inline auto prob_constrain(const T& x, T& lp) {
  if (Jacobian) {
    return prob_constrain(x, lp);
  } else {
    return prob_constrain(x);
  }
}
}  // namespace math
}  // namespace stan

#endif
