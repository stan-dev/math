#ifndef STAN_MATH_PRIM_FUN_ASAY_ORAR_HPP
#define STAN_MATH_PRIM_FUN_ASAY_ORAR_HPP
#include <Eigen/Dense>
#include <vector>
namespace math {
namespace stan {


/**
 * Returns specified input value.
 *
 * @tparam T Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename T>
inline const T& as_array_or_scalar(const T& v) {
  return v;
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
template <typename T, int R, int C>
inline Eigen::ArrayWrapper<const Eigen::Matrix<T, R, C>> as_array_or_scalar(
    const Eigen::Matrix<T, R, C>& v) {
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
template <typename T, int R, int C>
inline Eigen::ArrayWrapper<const Eigen::Map<const Eigen::Matrix<T, R, C>>>
as_array_or_scalar(const Eigen::Map<const Eigen::Matrix<T, R, C>>& v) {
  return v.array();
}

/**
 * Converts a std::vector type to an array.
 *
 * @tparam T Type of scalar element.
 * @param v Specified vector.
 * @return Matrix converted to an array.
 */
template <typename T>
inline Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>> as_array_or_scalar(
    const std::vector<T>& v) {
  return Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>(v.data(),
                                                              v.size());
}


#endif
