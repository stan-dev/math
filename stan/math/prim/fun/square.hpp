#ifndef STAN_MATH_PRIM_FUN_SQUARE_HPP
#define STAN_MATH_PRIM_FUN_SQUARE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the square of the specified argument.
 *
 * <p>\f$\mbox{square}(x) = x^2\f$.
 *
 * <p>The implementation of <code>square(x)</code> is just
 * <code>x * x</code>.  Given this, this method is mainly useful
 * in cases where <code>x</code> is not a simple primitive type,
 * particularly when it is an autodiff type.
 *
 * @param x Input to square.
 * @return Square of input.
 */
inline double square(double x) { return std::pow(x, 2); }

/**
 * Structure to wrap square() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return x squared.
 */
struct square_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return square(x);
  }
};

/**
 * Vectorized version of square().
 *
 * @tparam Container type of container
 * @param x container
 * @return Each value in x squared.
 */
template <
    typename Container, require_not_stan_scalar_t<Container>* = nullptr,
    require_not_container_st<std::is_arithmetic, Container>* = nullptr,
    require_not_var_matrix_t<Container>* = nullptr,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<Container>* = nullptr>
inline auto square(const Container& x) {
  return apply_scalar_unary<square_fun, Container>::apply(x);
}

/**
 * Version of square() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Each value in x squared.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto square(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().square(); });
}

}  // namespace math
}  // namespace stan

#endif
