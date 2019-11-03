#ifndef STAN_MATH_PRIM_MAT_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_MAT_FUN_DIVIDE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return specified matrix divided by specified scalar.
 * @tparam R Row type for matrix.
 * @tparam C Column type for matrix.
 * @param m Matrix.
 * @param c Scalar.
 * @return Matrix divided by scalar.
 */
template <int R, int C, typename T1, typename T2,
          typename = require_all_arithmetic_t<T1, T2>>
inline Eigen::Matrix<return_type_t<T1, T2>, R, C> divide(
    const Eigen::Matrix<T1, R, C>& m, T2 c) {
  return m / c;
}

}  // namespace math
}  // namespace stan
#endif
