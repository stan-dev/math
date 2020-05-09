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
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of the arguments' values
 */
template <typename T>
inline bool operator||(const var_value<T>& x, const var_value<T>& y) {
  return x.val() || y.val();
}

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam Arith An arithmetic type
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of first argument's value and second
 * argument
 */
template <typename T, typename Arith, require_arithmetic_t<Arith>...>
inline bool operator||(const var_value<T>& x, const Arith& y) {
  return x.val() || y;
}

/**
 * Return the logical disjunction of the values of the two
 * arguments as defined by <code>||</code>.
 *
 * @tparam Arith An arithmetic type
 * @param[in] x first argument
 * @param[in] y second argument
 * @return disjunction of first argument and the second
 * argument's value
 */
template <typename T, typename Arith, require_arithmetic_t<Arith>...>
inline bool operator||(const Arith& x, const var_value<T>& y) {
  return x || y.val();
}

}  // namespace math
}  // namespace stan
#endif
