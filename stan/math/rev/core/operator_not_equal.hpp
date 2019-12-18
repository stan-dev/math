#ifndef STAN_MATH_REV_CORE_OPERATOR_NOT_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_NOT_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Inequality operator comparing two variables' values (C++).
 *
   \f[
   \mbox{operator!=}(x, y) =
   \begin{cases}
     0 & \mbox{if } x = y\\
     1 & \mbox{if } x \neq y \\[6pt]
     0 & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam LHS value type of a var
 * @tparam RHS value type of a var
 * @param a First variable.
 * @param b Second variable.
 * @return True if the first variable's value is not the same as the
 * second's.
 */
template <typename LHS, typename RHS, require_all_var_t<LHS, RHS>...>
inline bool operator!=(LHS&& a, RHS&& b) {
  return a.val() != b.val();
}

/**
 * Inequality operator comparing a variable's value and a double
 * (C++).
 *
 * @tparam LHS value type of a var
 * @tparam RHS An arithmetic type
 * @param a First variable.
 * @param b Second value.
 * @return True if the first variable's value is not the same as the
 * second value.
 */
template <typename LHS, typename RHS,
require_var_t<LHS>..., require_arithmetic_t<RHS>...>
inline bool operator!=(LHS&& a, RHS b) { return a.val() != b; }

/**
 * Inequality operator comparing a double and a variable's value
 * (C++).
 *
 * @tparam LHS An arithmetic type
 * @tparam RHS value type of a var
 * @param a First value.
 * @param b Second variable.
 * @return True if the first value is not the same as the
 * second variable's value.
 */
template <typename LHS, typename RHS,
 require_arithmetic_t<LHS>..., require_var_t<RHS>...>
inline bool operator!=(LHS a, RHS&& b) { return a != b.val(); }

}  // namespace math
}  // namespace stan
#endif
