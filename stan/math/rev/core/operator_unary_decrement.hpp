#ifndef STAN_MATH_REV_CORE_OPERATOR_UNARY_DECREMENT_HPP
#define STAN_MATH_REV_CORE_OPERATOR_UNARY_DECREMENT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/callback_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Prefix decrement operator for variables (C++).
 *
 * Following C++, <code>(--a)</code> is defined to behave exactly as
 *
 * <code>a = a - 1.0)</code>
 *
 * does, but is faster and uses less memory.  In particular,
 * the result is an assignable lvalue.
 *
 * @param a Variable to decrement.
 * @return Reference the result of decrementing this input variable.
 */
inline var& operator--(var& a) {
  a = make_callback_var(a.val() - 1.0, [a](auto& vi) { a.adj() += vi.adj(); });
  return a;
}

/**
 * Postfix decrement operator for variables (C++).
 *
 * Following C++, the expression <code>(a--)</code> is defined to
 * behave like the sequence of operations
 *
 * <code>var temp = a;  a = a - 1.0;  return temp;</code>
 *
 * @param a Variable to decrement.
 * @return Input variable.
 */
inline var operator--(var& a, int /*dummy*/) {
  var temp(a);
  a = make_callback_var(a.val() - 1.0, [a](auto& vi) { a.adj() += vi.adj(); });
  return temp;
}

}  // namespace math
}  // namespace stan
#endif
