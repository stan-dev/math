#ifndef STAN_MATH_REV_CORE_OPERATOR_LESS_THAN_HPP
#define STAN_MATH_REV_CORE_OPERATOR_LESS_THAN_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Less than operator comparing variables' values (C++).
 *
 \f[
   \mbox{operator\textless}(x, y) =
   \begin{cases}
     0 & \mbox{if } x \geq y \\
     1 & \mbox{if } x < y \\[6pt]
     0 & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
 \f]
 * @tparam Var1 value type of a var
 * @tparam Var2 value type of a var
 * @param a First variable.
 * @param b Second variable.
 * @return True if first variable's value is less than second's.
 */
template <typename Var1, typename Var2, require_all_var_t<Var1, Var2>...>
inline bool operator<(Var1&& a, Var2&& b) { return a.val() < b.val(); }

/**
 * Less than operator comparing variable's value and a double
 * (C++).
 *
 * @tparam Arith An arithmetic type
 * @tparam Var value type of a var
 * @param a First variable.
 * @param b Second value.
 * @return True if first variable's value is less than second value.
 */
 template <typename Var, typename Arith,
  require_var_t<Var>..., require_arithmetic_t<Arith>...>
inline bool operator<(Var&& a, Arith b) { return a.val() < b; }

/**
 * Less than operator comparing a double and variable's value
 * (C++).
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param a First value.
 * @param b Second variable.
 * @return True if first value is less than second variable's value.
 */
 template <typename Var, typename Arith,
  require_var_t<Var>..., require_arithmetic_t<Arith>...>
inline bool operator<(Arith a, Var&& b) { return a < b.val(); }

}  // namespace math
}  // namespace stan
#endif
