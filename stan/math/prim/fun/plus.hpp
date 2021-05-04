#ifndef STAN_MATH_PRIM_FUN_PLUS_HPP
#define STAN_MATH_PRIM_FUN_PLUS_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the unary plus of the input.
 *
 * @tparam T Type of input.
 * @param x input.
 * @return result of unary plus of the input.
 */
template <typename T>
inline T plus(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan

#endif
