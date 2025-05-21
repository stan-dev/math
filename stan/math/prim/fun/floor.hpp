#ifndef STAN_MATH_PRIM_FUN_FLOOR_HPP
#define STAN_MATH_PRIM_FUN_FLOOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto floor(const T x) {
  return std::floor(x);
}

template <typename T, require_complex_bt<std::is_arithmetic, T>* = nullptr>
inline auto floor(const T x) {
  return std::floor(x);
}
/**
 * Structure to wrap `floor()` so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Greatest integer <= x.
 */
struct floor_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return floor(x);
  }
};

/**
 * Returns the elementwise `floor()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return Greatest integer <= each value in x.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr>
inline auto floor(const Container& x) {
  return apply_scalar_unary<floor_fun, Container>::apply(x);
}

/**
 * Version of `floor()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Greatest integer <= each value in x.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr>
inline auto floor(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().floor(); });
}

}  // namespace math
}  // namespace stan

#endif
