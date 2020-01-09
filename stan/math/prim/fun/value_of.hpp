#ifndef STAN_MATH_PRIM_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_FUN_VALUE_OF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the value of the specified scalar argument
 * converted to a double value.
 *
 * <p>See the <code>primitive_value</code> function to
 * extract values without casting to <code>double</code>.
 *
 * <p>This function is meant to cover the primitive types. For
 * types requiring pass-by-reference, this template function
 * should be specialized.
 *
 * @tparam T type of scalar.
 * @param x scalar to convert to double
 * @return value of scalar cast to double
 */
template <typename T>
inline double value_of(const T x) {
  return static_cast<double>(x);
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x value
 * @return input value
 */
template <>
inline double value_of<double>(double x) {
  return x;
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x value
 * @return input value
 */
inline int value_of(int x) { return x; }

/**
 * Convert a std::vector of type T to a std::vector of
 * child_type<T>::type.
 *
 * @tparam T Scalar type in std::vector
 * @param[in] x std::vector to be converted
 * @return std::vector of values
 **/
template <typename T>
inline std::vector<typename child_type<T>::type> value_of(
    const std::vector<T>& x) {
  size_t x_size = x.size();
  std::vector<typename child_type<T>::type> result(x_size);
  for (size_t i = 0; i < x_size; i++) {
    result[i] = value_of(x[i]);
  }
  return result;
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified std::vector.
 * @return Specified std::vector.
 */
inline const std::vector<double>& value_of(const std::vector<double>& x) {
  return x;
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified std::vector.
 * @return Specified std::vector.
 */
inline const std::vector<int>& value_of(const std::vector<int>& x) { return x; }

/**
 * Convert a matrix of type T to a matrix of doubles.
 *
 * T must implement value_of. See
 * test/math/fwd/fun/value_of.cpp for fvar and var usage.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows in the matrix, can be Eigen::Dynamic
 * @tparam C number of columns in the matrix, can be Eigen::Dynamic
 *
 * @param[in] M Matrix to be converted
 * @return Matrix of values
 **/
template <typename T, int R, int C>
inline Eigen::Matrix<typename child_type<T>::type, R, C> value_of(
    const Eigen::Matrix<T, R, C>& M) {
  Eigen::Matrix<typename child_type<T>::type, R, C> Md(M.rows(), M.cols());
  for (int j = 0; j < M.cols(); j++) {
    for (int i = 0; i < M.rows(); i++) {
      Md(i, j) = value_of(M(i, j));
    }
  }
  return Md;
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @tparam R number of rows in the matrix, can be Eigen::Dynamic
 * @tparam C number of columns in the matrix, can be Eigen::Dynamic
 *
 * @param x Specified matrix.
 * @return Specified matrix.
 */
template <int R, int C>
inline const Eigen::Matrix<double, R, C>& value_of(
    const Eigen::Matrix<double, R, C>& x) {
  return x;
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @tparam R number of rows in the matrix, can be Eigen::Dynamic
 * @tparam C number of columns in the matrix, can be Eigen::Dynamic
 *
 * @param x Specified matrix.
 * @return Specified matrix.
 */
template <int R, int C>
inline const Eigen::Matrix<int, R, C>& value_of(
    const Eigen::Matrix<int, R, C>& x) {
  return x;
}

}  // namespace math
}  // namespace stan

#endif
