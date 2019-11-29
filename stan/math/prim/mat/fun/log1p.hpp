#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap log1p() so it can be vectorized.
 */
struct log1p_fun {
  /**
   * Return the inverse hypberbolic cosine of the specified argument.
   *
   * @param x Argument.
   * @return Inverse hyperbolic cosine of the argument.
   * @tparam T Argument type.
   */
  template <typename T>
  static inline T fun(const T& x) {
    return log1p(x);
  }
};

/**
 * Return the elementwise application of <code>log1p()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T Container type.
 * @param x Container.
 * @return Elementwise log1p of members of container.
 */
template <
    typename T /*, typename = require_not_eigen_vt<std::is_arithmetic, T>*/>
inline auto log1p(const T& x) {
  return apply_scalar_unary<log1p_fun, T>::apply(x);
}

// TODO(Tadej): This is inconsistent with stan's scalar version (which throws on
// value <-1).
///**
// * Version of log1p() that accepts Eigen Matrix ar matrix expressions.
// * @tparam Derived derived type of x
// * @param x Matrix or matrix expression
// * @return Elementwise log1p of members of container.
// */
// template <typename Derived, typename = require_eigen_vt<std::is_arithmetic,
// Derived>>  inline auto log1p(const Eigen::MatrixBase<Derived>& x){
//  return x.derived().array().log1p().matrix();
//}

}  // namespace math
}  // namespace stan

#endif
