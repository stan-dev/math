#ifndef STAN_MATH_PRIM_FUN_LMULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_LMULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the first argument times the log of the second argument. The
 * result is 0 if both arguments are 0.  The funcgtion is defined by
 * `lmultiply(x, y) = x * log(y)` if `x` or `y` is non-zero and
 * `lmultiply(0, 0) = 0` otherwise.
 *
 * @tparam T1 type of the first argument
 * @tparam T2 type of the second argument
 * @param a first argument
 * @param b second argument
 * @return the first argument times the log of the second argument
 */
template <typename T1, typename T2>
inline auto lmultiply(T1&& a, T2&& b) {
  if constexpr (is_kernel_expression<T1>::value
                || is_kernel_expression<T2>::value) {
    return multiply_log(std::forward<T1>(a), std::forward<T2>(b));
  } else {
    return make_holder(
        [](auto&& a, auto&& b) {
          return multiply_log(std::forward<decltype(a)>(a),
                              std::forward<decltype(b)>(b));
        },
        std::forward<T1>(a), std::forward<T2>(b));
  }
}
}  // namespace math
}  // namespace stan

#endif
