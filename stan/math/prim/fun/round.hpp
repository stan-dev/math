#ifndef STAN_MATH_PRIM_FUN_ROUND_HPP
#define STAN_MATH_PRIM_FUN_ROUND_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap `round()` so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x argument variable
 * @return Rounded value of x.
 */
struct round_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    using std::round;
    return round(x);
  }
};

/**
 * Vectorized version of `round()`.
 *
 * @tparam Container type of container
 * @param x container
 * @return Rounded value of each value in x.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto round(const Container& x) {
  return apply_scalar_unary<round_fun, Container>::apply(x);
}

/**
 * Version of `round()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Rounded value of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto round(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().round(); });
}

}  // namespace math
}  // namespace stan

#endif
