#ifndef STAN_MATH_PRIM_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_FUN_INV_SQRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_stan_scalar_t<T>* = nullptr>
inline auto inv_sqrt(T x) {
  using std::sqrt;
  return inv(sqrt(x));
}
/**
 * Structure to wrap `1 / sqrt(x)` so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return inverse square root of x.
 */
struct inv_sqrt_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return inv_sqrt(x);
  }
};

/**
 * Return the elementwise `1 / sqrt(x)}` of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return inverse square root of each value in x.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr,
          require_not_stan_scalar_t<Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto inv_sqrt(const Container& x) {
  return apply_scalar_unary<inv_sqrt_fun, Container>::apply(x);
}

/**
 * Version of `inv_sqrt()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return inverse square root each variable in the container.
 */
template <typename Container, require_not_var_matrix_t<Container>* = nullptr,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto inv_sqrt(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().rsqrt(); });
}

}  // namespace math
}  // namespace stan

#endif
