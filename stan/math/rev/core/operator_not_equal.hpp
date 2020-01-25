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
 * @param a First variable.
 * @param b Second variable.
 * @return True if the first variable's value is not the same as the
 * second's.
 */
inline bool operator!=(var a, var b) { return a.val() != b.val(); }

/**
 * Inequality operator comparing a variable's value and a double
 * (C++).
 *
 * @tparam Arith An arithmetic type
 * @param a First variable.
 * @param b Second value.
 * @return True if the first variable's value is not the same as the
 * second value.
 */
template <typename Arith, require_arithmetic_t<Arith>...>
inline bool operator!=(var a, Arith b) {
  return a.val() != b;
}

/**
 * Inequality operator comparing a double and a variable's value
 * (C++).
 *
 * @tparam Arith An arithmetic type
 * @param a First value.
 * @param b Second variable.
 * @return True if the first value is not the same as the
 * second variable's value.
 */
template <typename Arith, require_arithmetic_t<Arith>...>
inline bool operator!=(Arith a, var b) {
  return a != b.val();
}

}  // namespace math
}  // namespace stan
#endif
