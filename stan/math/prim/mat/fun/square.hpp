#ifndef STAN_MATH_PRIM_MAT_FUN_SQUARE_HPP
#define STAN_MATH_PRIM_MAT_FUN_SQUARE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

namespace stan {
namespace math {

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
 * @tparam T type of container
 * @param x container
 * @return Each value in x squared.
 */
template <typename T>
inline auto square(const T& x) {
  return apply_scalar_unary<square_fun, T>::apply(x);
}

/**
 * Version of square() that accepts Eigen Matrix or matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Each value in x squared.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto square(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().square().matrix();
}

}  // namespace math
}  // namespace stan

#endif
