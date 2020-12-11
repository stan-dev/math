#ifndef STAN_MATH_PRIM_FUN_MINUS_HPP
#define STAN_MATH_PRIM_FUN_MINUS_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the negation of the specified scalar or matrix.
 *
 * @tparam T Type of subtrahend.
 * @param x Subtrahend.
 * @return Negation of subtrahend.
 */
template <typename T>
inline auto minus(const T& x) {
  return -x;
}

}  // namespace math
}  // namespace stan

#endif
