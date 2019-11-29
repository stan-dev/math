#ifndef STAN_MATH_PRIM_MAT_FUN_INV_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return 1 / x.
 */
struct inv_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv(x);
  }
};

/**
 * Vectorized version of inv().
 * @param x Container.
 * @tparam T Container type.
 * @return 1 divided by each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto inv(const T& x) {
  return apply_scalar_unary<inv_fun, T>::apply(x);
}

/**
 * Version of inv() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return 1 divided by each value in x.
 */
template <typename Derived, typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv(const Eigen::MatrixBase<Derived>& x){
  return x.derived().array().inverse().matrix();
}

/**
 * Version of inv() that accepts Eigen Array ar array expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return 1 divided by each value in x.
 */
template <typename Derived, typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv(const Eigen::ArrayBase<Derived>& x){
  return x.derived().inverse();
}

}  // namespace math
}  // namespace stan

#endif
