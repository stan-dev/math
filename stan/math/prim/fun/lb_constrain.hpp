#ifndef STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of scalar
 * @tparam L type of lower bound
 * @param[in] x Unconstrained scalar input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, L> lb_constrain(const T& x, const L& lb) {
  using std::exp;
  if (unlikely(is_negative_infinity(value_of(lb)))) {
    return identity_constrain(x);
  }
  return exp(x) + lb;
}

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of Matrix
 * @tparam L type of lower bound
 * @param[in] x Unconstrained Matrix input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_matrix_t<T>* = nullptr>
inline auto lb_constrain(T&& x, const L& lb) {
  if (unlikely(is_negative_infinity(value_of(lb)))) {
    return identity_constrain(std::forward<T>(x)).eval();
  }
  return add(exp(x), lb).eval();
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * If the lower bound is negative infinity, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam T type of scalar
 * @tparam L type of lower bound
 * @tparam S type of log probability
 * @param[in] x unconstrained scalar input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, typename S, require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, L> lb_constrain(const T& x, const L& lb, S& lp) {
  using std::exp;
  if (unlikely(is_negative_infinity(value_of(lb)))) {
    return identity_constrain(x, lp);
  }
  lp += x;
  return exp(x) + lb;
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * If the lower bound is negative infinity, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam T Type of Matrix
 * @tparam L type of lower bound
 * @tparam S type of log probability
 * @param[in] x unconstrained Matrix input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, typename S, require_matrix_t<T>* = nullptr>
inline auto lb_constrain(T&& x, const L& lb, S& lp) {
  if (unlikely(is_negative_infinity(value_of(lb)))) {
    return identity_constrain(std::forward<T>(x), lp).eval();
  }
  lp += sum(x);
  return add(exp(x), lb).eval();
}

}  // namespace math

}  // namespace stan

#endif
