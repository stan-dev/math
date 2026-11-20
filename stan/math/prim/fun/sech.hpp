#ifndef STAN_MATH_PRIM_FUN_SECH_HPP
#define STAN_MATH_PRIM_FUN_SECH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/cosh.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the hyperbolic secant `1 / cosh(x)` of the argument, elementwise
 * for containers. For large `|x|`, `cosh(x)` overflows to infinity and the
 * result is 0.
 *
 * @tparam T scalar, container, or `var_value` with inner Eigen type
 * @param x argument
 * @return hyperbolic secant of the argument
 */
template <typename T>
inline auto sech(T&& x) {
  return inv(cosh(std::forward<T>(x)));
}

}  // namespace math
}  // namespace stan

#endif
