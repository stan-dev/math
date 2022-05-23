#ifndef STAN_MATH_PRIM_FUN_INV_HPP
#define STAN_MATH_PRIM_FUN_INV_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap 1.0 / x so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 / x.
 */
struct inv_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return 1.0 / x;
  }
};

/**
 * Return the elementwise 1.0 / x of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return 1 divided by each value in x.
 */
template <
    typename T, require_not_container_st<std::is_arithmetic, T>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto inv(const T& x) {
  return apply_scalar_unary<inv_fun, T>::apply(x);
}

/**
 * Version of \c inv() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return 1 divided by each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto inv(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().inverse(); });
}

}  // namespace math
}  // namespace stan

#endif
