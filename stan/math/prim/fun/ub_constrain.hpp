#ifndef STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_var_value_if.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

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
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @return matrix constrained to have upper bound
 */
template <typename T, typename L>
inline auto ub_constrain(const T& x, const L& ub) {
  auto&& x_ref = to_ref(x);
  auto&& ub_ref = to_ref(ub);
  if (is_positive_infinity(ub_ref)) {
    return identity_constrain(x, ub_ref);
  } else {
    return eval(subtract(ub_ref, exp(x_ref)));
  }
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
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename T, typename L>
inline auto ub_constrain(const T& x, const L& ub, return_type_t<T, L>& lp) {
  auto&& x_ref = to_ref(x);
  auto&& ub_ref = to_ref(ub);
  if (is_positive_infinity(ub_ref)) {
    return identity_constrain(x, ub_ref);
  } else {
    lp += sum(x_ref);
    return eval(subtract(ub_ref, exp(x_ref)));
  }
}

}  // namespace math
}  // namespace stan

#endif
