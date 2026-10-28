#ifndef STAN_MATH_PRIM_FUN_MULTIPLY_LOG_HPP
#define STAN_MATH_PRIM_FUN_MULTIPLY_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/lmultiply.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the first argument times the log of the second argument.
 * This is the deprecated name of `lmultiply`, kept because code generated
 * by released stanc3 versions calls it; it can be removed only in a major
 * release. Library code calls `lmultiply`, which has the definition and
 * the derivatives.
 *
 * @tparam T1 type of the first argument
 * @tparam T2 type of the second argument
 * @param a first argument
 * @param b second argument
 * @return `lmultiply(a, b)`
 */
template <typename T1, typename T2>
inline auto multiply_log(T1&& a, T2&& b) {
  return lmultiply(std::forward<T1>(a), std::forward<T2>(b));
}

}  // namespace math
}  // namespace stan

#endif
