#ifndef STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
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
 * @tparam T A scalar or matrix.
 * @tparam L A scalar or matrix.
 * @param[in] x Unconstrained input
 * @param[in] lb Lower bound
 * @return Constrained matrix
 */
template <typename T, typename L>
inline auto lb_constrain(const T& x, const L& lb) {
  auto&& x_ref = to_ref(x);
  auto&& lb_ref = to_ref(lb);
  if (is_negative_infinity(lb_ref)) {
    return identity_constrain(x_ref, lb_ref);
  } else {
    return eval(add(exp(x_ref), lb_ref));
  }
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * @tparam T A scalar or matrix.
 * @tparam L A scalar or matrix.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  auto&& x_ref = to_ref(x);
  auto&& lb_ref = to_ref(lb);
  if (is_negative_infinity(lb_ref)) {
    return identity_constrain(x_ref, lb_ref);
  } else {
    lp += sum(x_ref);
    return eval(add(exp(x_ref), lb_ref));
  }
}

}  // namespace math
}  // namespace stan

#endif
