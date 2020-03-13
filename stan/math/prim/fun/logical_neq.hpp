#ifndef STAN_MATH_PRIM_FUN_LOGICAL_NEQ_HPP
#define STAN_MATH_PRIM_FUN_LOGICAL_NEQ_HPP

#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * Return 1 if the first argument is unequal to the second.
 * Equivalent to <code>x1 != x2</code>.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param x1 first argument
 * @param x2 second argument
 * @return <code>true</code> iff <code>x1 != x2</code>
 */
template <typename T1, typename T2>
inline int logical_neq(const T1 x1, const T2 x2) {
  return x1 != x2;
}

}  // namespace math
}  // namespace stan
#endif
