#ifndef STAN_MATH_PRIM_MAT_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_MAT_FUN_DIVIDE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return matrix divided by scalar.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @tparam T1 type of elements in the matrix
 * @tparam T2 type of scalar
 * @param[in] m specified matrix
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
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
