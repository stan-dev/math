#ifndef STAN_MATH_REV_CORE_OPERATOR_GREATER_THAN_HPP
#define STAN_MATH_REV_CORE_OPERATOR_GREATER_THAN_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Greater than operator comparing variables' values (C++).
 *
   \f[
   \mbox{operator\textgreater}(x, y) =
   \begin{cases}
     0 & \mbox{if } x \leq y\\
     1 & \mbox{if } x > y \\[6pt]
     0 & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam LHS value type of a var
 * @tparam RHS value type of a var
 * @param a First variable.
 * @param b Second variable.
 * @return True if first variable's value is greater than second's.
 */
template <typename LHS, typename RHS, require_all_var_t<LHS, RHS>...>
inline bool operator>(LHS&& a, RHS&& b) { return a.val() > b.val(); }

/**
 * Greater than operator comparing variable's value and double
 * (C++).
 *
 * @tparam LHS value type of a var
 * @tparam RHS An arithmetic type
 * @param a First variable.
 * @param b Second value.
 * @return True if first variable's value is greater than second value.
 */
 template <typename LHS, typename RHS,
  require_var_t<LHS>..., require_arithmetic_t<RHS>...>
inline bool operator>(LHS&& a, RHS b) { return a.val() > b; }

/**
 * Greater than operator comparing a double and a variable's value
 * (C++).
 *
 * @tparam LHS an Arithmetic type 
 * @tparam RHS value type of a var
 * @param a First value.
 * @param b Second variable.
 * @return True if first value is greater than second variable's value.
 */
template <typename LHS, typename RHS,
 require_arithmetic_t<LHS>..., require_var_t<RHS>...>
inline bool operator>(LHS a, RHS&& b) { return a > b.val(); }

}  // namespace math
}  // namespace stan
#endif
