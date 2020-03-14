#ifndef STAN_MATH_PRIM_FUN_EXP_HPP
#define STAN_MATH_PRIM_FUN_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {
/**
 * Return the natural exponential of the specified argument.  This
 * version is required to disambiguate <code>exp(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural exponential of argument.
 */
inline double exp(int x) { return std::exp(x); }

/**
 * Structure to wrap <code>exp()</code> so that it can be
 * vectorized.
 */
struct exp_fun {
  /**
   * Return the exponential of the specified scalar argument.
   *
   * @tparam T type of argument
   * @param[in] x argument
   * @return Exponential of argument.
   */
  template <typename T>
  static inline T fun(const T& x) {
    using std::exp;
    return exp(x);
  }
};

/**
 * Return the elementwise exponentiation of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam Container type of container
 * @param[in] x container
 * @return Elementwise application of exponentiation to the argument.
 */
template <
    typename Container,
    require_not_container_st<is_container, std::is_arithmetic, Container>...>
inline auto exp(const Container& x) {
  return apply_scalar_unary<exp_fun, Container>::apply(x);
}

/**
 * Version of exp() that accepts std::vectors, Eigen Matrix/Array objects
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
      x, [](auto&& v) { return v.array().exp(); });
}

}  // namespace math
}  // namespace stan

#endif
