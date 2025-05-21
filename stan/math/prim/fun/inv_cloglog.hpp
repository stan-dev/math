#ifndef STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_PRIM_FUN_INV_CLOGLOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse complementary log-log function.
 *
 * The function is defined by
 *
 * <code>inv_cloglog(x) = 1 - exp(-exp(x))</code>.
 *
 * This function can be used to implement the inverse link
 * function for complementary-log-log regression.
 *
 *
   \f[
   \mbox{inv\_cloglog}(y) =
   \begin{cases}
     \mbox{cloglog}^{-1}(y) & \mbox{if } -\infty\leq y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } y = \textrm{NaN}
   \end{cases}
   \f]
   \f[
   \frac{\partial\, \mbox{inv\_cloglog}(y)}{\partial y} =
   \begin{cases}
     \frac{\partial\, \mbox{cloglog}^{-1}(y)}{\partial y} & \mbox{if }
 -\infty\leq y\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } y = \textrm{NaN}
   \end{cases}
   \f]
   \f[
   \mbox{cloglog}^{-1}(y) = 1 - \exp \left( - \exp(y) \right)
   \f]
   \f[
   \frac{\partial \, \mbox{cloglog}^{-1}(y)}{\partial y} = \exp(y-\exp(y))
   \f]
 *
 * @param x Argument.
 * @return Inverse complementary log-log of the argument.
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto inv_cloglog(const T x) {
  return 1. - std::exp(-std::exp(x));
}

/**
 * The inverse complementary log-log function.
 *
 * The function is defined by
 *
 * <code>inv_cloglog(x) = 1 - exp(-exp(x))</code>.
 *
 * This function can be used to implement the inverse link
 * function for complementary-log-log regression.
 *
 * @param x Argument.
 * @return Inverse complementary log-log of the argument.
 */
template <typename T, require_complex_t<T>* = nullptr>
inline auto inv_cloglog(const T& x) {
  return 1. - exp(-exp(x));
}

/**
 * Structure to wrap inv_cloglog() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 - exp(-exp(x)).
 */
struct inv_cloglog_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return inv_cloglog(x);
  }
};

/**
 * Vectorized version of inv_cloglog().
 *
 * @tparam Container type of container
 * @param x container
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr>
inline auto inv_cloglog(const Container& x) {
  return apply_scalar_unary<inv_cloglog_fun, Container>::apply(x);
}

/**
 * Version of inv_cloglog() that accepts std::vectors, Eigen Matrix/Array
 * objects or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto inv_cloglog(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return 1 - (-v.array().exp()).exp(); });
}

}  // namespace math
}  // namespace stan

#endif
