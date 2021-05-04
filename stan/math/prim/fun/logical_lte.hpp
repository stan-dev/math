#ifndef STAN_MATH_PRIM_FUN_LOGICAL_LTE_HPP
#define STAN_MATH_PRIM_FUN_LOGICAL_LTE_HPP

namespace stan {
namespace math {

/**
 * Return 1 if the first argument is less than or equal to the second.
 * Equivalent to <code>x1 &lt;= x2</code>.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param x1 first argument
 * @param x2 second argument
 * @return <code>true</code> iff <code>x1 &lt;= x2</code>
 */
template <typename T1, typename T2>
inline bool logical_lte(const T1 x1, const T2 x2) {
  return x1 <= x2;
}

}  // namespace math
}  // namespace stan

#endif
