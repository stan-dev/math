#ifndef STAN_MATH_REV_CORE_OPERATOR_LOGICAL_AND_HPP
#define STAN_MATH_REV_CORE_OPERATOR_LOGICAL_AND_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the logical conjunction of the values of the two
 * arguments as defined by <code>&amp;&amp;</code>.
 *
 * @tparam LHS value type of a var
 * @tparam RHS value type of a var
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjuntion of the argument's values
 */
template <typename LHS, typename RHS, require_all_var_t<LHS, RHS>...>
inline bool operator&&(LHS&& x, RHS&& y) {
  return x.val() && y.val();
}

/**
 * Return the logical conjunction of the values of the two
 * arguments as defined by <code>&amp;&amp;</code>.
 *
 * @tparam RHS An arithmetic type
 * @tparam LHS value type of a var
 * @param[in] x first argument
 * @param[in] y second argument
 * @return conjunction of first argument's value and second
 * argument
 */
template <typename LHS, typename RHS,
 require_var_t<LHS>..., require_arithmetic_t<RHS>...>
inline bool operator&&(LHS&& x, RHS y) {
  return x.val() && y;
}

/**
 * Return the logical conjunction of the values of the two
 * arguments as defined by <code>&amp;&amp;</code>.
 *
 * @tparam LHS An arithmetic type
 * @tparam RHS value type of a var
 * @param[in] x first argument
 * @param[in] y second argument
 * @return conjunction of first argument and second argument's
 * value
 */
template <typename LHS, typename RHS,
 require_arithmetic_t<LHS>..., require_var_t<RHS>...>
inline bool operator&&(LHS x, RHS&& y) {
  return x && y.val();
}

}  // namespace math
}  // namespace stan
#endif
