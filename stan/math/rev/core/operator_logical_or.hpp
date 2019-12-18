#ifndef STAN_MATH_REV_CORE_OPERATOR_LOGICAL_OR_HPP
#define STAN_MATH_REV_CORE_OPERATOR_LOGICAL_OR_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam Var1 value type of a var
 * @tparam Var2 value type of a var
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjuntion of the argument's values
 */
template <typename Var1, typename Var2, require_all_var_t<Var1, Var2>...>
inline bool operator||(Var1&& x, Var2&& y) {
  return x.val() || y.val();
}

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of first argument's value and second
 * argument
 */
 template <typename Var, typename Arith,
  require_var_t<Var>..., require_arithmetic_t<Arith>...>
inline bool operator||(Var&& x, Arith y) {
  return x.val() || y;
}

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of first argument and the second
 * argument's value
 */
 template <typename Var, typename Arith,
  require_var_t<Var>..., require_arithmetic_t<Arith>...>
inline bool operator||(Arith x, Var&& y) {
  return x || y.val();
}

}  // namespace math
}  // namespace stan
#endif
