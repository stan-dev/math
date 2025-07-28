#ifndef STAN_MATH_PRIM_FUN_FABS_HPP
#define STAN_MATH_PRIM_FUN_FABS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto fabs(T&& x) {
  return std::abs(x);
}

template <typename T, require_complex_bt<std::is_arithmetic, T>* = nullptr>
inline auto fabs(T&& x) {
  return std::hypot(x.real(), x.imag());
}

/**
 * Structure to wrap `fabs()` so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Absolute value of x.
 */
struct fabs_fun {
  template <typename T>
  static inline auto fun(T&& x) {
    return fabs(std::forward<T>(x));
  }
};

/**
 * Returns the elementwise `fabs()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return Absolute value of each value in x.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr,
          require_not_stan_scalar_t<Container>* = nullptr,
          require_container_t<Container>* = nullptr>
inline auto fabs(Container&& x) {
  return apply_scalar_unary<fabs_fun, Container>::apply(
      std::forward<Container>(x));
}

/**
 * Version of `fabs()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Absolute value of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto fabs(Container&& x) {
  return apply_vector_unary<Container>::apply(
      std::forward<Container>(x), [](auto&& v) { return v.array().abs(); });
}

}  // namespace math
}  // namespace stan

#endif
