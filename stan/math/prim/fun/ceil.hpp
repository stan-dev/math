#ifndef STAN_MATH_PRIM_FUN_CEIL_HPP
#define STAN_MATH_PRIM_FUN_CEIL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto ceil(T&& x) {
  return std::ceil(x);
}

/**
 * Structure to wrap `ceil()` so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Least integer >= x.
 */
struct ceil_fun {
  template <typename T>
  static inline auto fun(T&& x) {
    return ceil(std::forward<T>(x));
  }
};

/**
 * Returns the elementwise `ceil()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return Least integer >= each value in x.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr,
          require_container_t<Container>* = nullptr>
inline auto ceil(Container&& x) {
  return apply_scalar_unary<ceil_fun, Container>::apply(
      std::forward<Container>(x));
}

/**
 * Version of `ceil()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Least integer >= each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr>
inline auto ceil(Container&& x) {
  return apply_vector_unary<Container>::apply(
      std::forward<Container>(x), [](auto&& v) { return v.array().ceil(); });
}

}  // namespace math
}  // namespace stan

#endif
