#ifndef STAN_MATH_PRIM_FUN_LOG10_HPP
#define STAN_MATH_PRIM_FUN_LOG10_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap log10() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Log base-10 of x.
 */
struct log10_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    using std::log10;
    return log10(x);
  }
};

/**
 * Vectorized version of log10().
 *
 * @tparam Container type of container
 * @param x container
 * @return Log base-10 applied to each value in x.
 */
template <
    typename Container, require_not_var_matrix_t<Container>* = nullptr,
    require_not_container_st<std::is_arithmetic, Container>* = nullptr,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<Container>* = nullptr>
inline auto log10(const Container& x) {
  return apply_scalar_unary<log10_fun, Container>::apply(x);
}

/**
 * Version of log10() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Log base-10 of each variable in the container.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto log10(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().log10(); });
}

namespace internal {
/**
 * Return the base 10 logarithm of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return base 10 logarithm of the argument
 */
template <typename V>
inline std::complex<V> complex_log10(const std::complex<V>& z) {
  static constexpr double inv_log_10 = 1.0f / LOG_TEN;
  return log(z) * inv_log_10;
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
