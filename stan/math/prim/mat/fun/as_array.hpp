#ifndef STAN_MATH_PRIM_MAT_FUN_AS_ARRAY_HPP
#define STAN_MATH_PRIM_MAT_FUN_AS_ARRAY_HPP

#include <Eigen\Dense>

namespace stan {
namespace math {

/**
 * Converts a matrix type to an array.
 *
 * @tparam T Type of scalar element.
 * @tparam R Row type of input matrix.
 * @tparam C Column type of input matrix.
 * @param v Specified matrix.
 * @return Matrix converted to an array.
 */
template<typename T, int R, int C>
inline Eigen::ArrayWrapper<Eigen::Matrix<T, R, C>> as_array(Eigen::Matrix<T, R, C>& v) {
  return v.array();
}

/**
 * Converts a matrix type to an array.
 *
 * @tparam T Type of scalar element.
 * @tparam R Row type of input matrix.
 * @tparam C Column type of input matrix.
 * @param v Specified matrix.
 * @return Matrix converted to an array.
 */
template<typename T, int R, int C>
inline Eigen::ArrayWrapper<const Eigen::Matrix<T, R, C>> as_array(const Eigen::Matrix<T, R, C>& v) {
  return v.array();
}

}
}

#endif
