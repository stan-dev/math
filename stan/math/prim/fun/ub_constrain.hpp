#ifndef STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @param[in] x free scalar.
 * @param[in] ub upper bound
 * @return scalar constrained to have upper bound
 */
template <typename T, typename U, require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, U> ub_constrain(const T& x, const U& ub) {
  using std::exp;
  if (unlikely(is_infinity(value_of(ub)))) {
    return identity_constrain(x);
  }
  return ub - exp(x);
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @return matrix constrained to have upper bound
 */
template <typename T, typename U, require_matrix_t<T>* = nullptr>
inline auto ub_constrain(T&& x, const U& ub) {
  if (unlikely(is_infinity(value_of(ub)))) {
    return identity_constrain(std::forward<T>(x));
  }
  return subtract(ub, exp(x)).eval();
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, lp)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename T, typename U, typename S,
          require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, U> ub_constrain(const T& x, const U& ub, S& lp) {
  using std::exp;
  if (unlikely(is_infinity(value_of(ub)))) {
    return identity_constrain(x, lp);
  }
  lp += x;
  return ub - exp(x);
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, lp)</code>.
 *
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free Matrix
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return matrix constrained to have upper bound
 */
template <typename T, typename U, typename S, require_matrix_t<T>* = nullptr>
inline auto ub_constrain(T&& x, const U& ub, S& lp) {
  if (unlikely(is_infinity(value_of(ub)))) {
    return identity_constrain(std::forward<T>(x), lp);
  }
  lp += sum(x);
  return subtract(ub, exp(x)).eval();
}

}  // namespace math

}  // namespace stan

#endif
