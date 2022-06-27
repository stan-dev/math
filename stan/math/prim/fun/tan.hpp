#ifndef STAN_MATH_PRIM_FUN_TAN_HPP
#define STAN_MATH_PRIM_FUN_TAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/i_times.hpp>
#include <stan/math/prim/fun/tanh.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap `tan()` so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Tangent of x.
 */
struct tan_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    using std::tan;
    return tan(x);
  }
};

/**
 * Vectorized version of `tan()`.
 *
 * @tparam Container type of container
 * @param x angles in radians
 * @return Tangent of each value in x.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto tan(const Container& x) {
  return apply_scalar_unary<tan_fun, Container>::apply(x);
}

/**
 * Version of `tan()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Tangent of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto tan(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().tan(); });
}

namespace internal {
/**
 * Return the tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return tangent of the argument
 */
template <typename V>
inline std::complex<V> complex_tan(const std::complex<V>& z) {
  return neg_i_times(tanh(i_times(z)));
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
