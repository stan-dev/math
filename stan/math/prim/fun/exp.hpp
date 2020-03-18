#ifndef STAN_MATH_PRIM_FUN_EXP_HPP
#define STAN_MATH_PRIM_FUN_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the exponential of the specified scalar argument.
 *
 * @tparam T type of argument
 * @param[in] x argument
 * @return Exponential of argument.
 */
template <typename T, require_arithmetic_t<T>...>
auto exp(const T& x) {
  using std::exp;
  return exp(x);
}

/**
 * Version of `exp()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Elementwise application of exponentiation to the argument.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto exp(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().exp(); });
}

}  // namespace math
}  // namespace stan

#endif
