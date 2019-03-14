#ifndef STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP

#include <vector>
#include <cstddef>

#include <stan/math/prim/fun/Eigen.hpp>

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
 * @tparam T Type of scalar.
 * @param x Scalar to convert to double.
 * @return Value of scalar cast to a double.
 */
template <typename T>
inline double value_of_rec(const T x) {
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
 * @param x Specified value.
 * @return Specified value.
 */
template <>
inline double value_of_rec<double>(double x) {
  return x;
}

/**
 * Convert a std::vector of type T to a std::vector of doubles.
 *
 * T must implement value_of_rec. See
 * test/math/fwd/mat/fun/value_of_rec.cpp for fvar and var usage.
 *
 * @tparam T Scalar type in std::vector
 * @param[in] x std::vector to be converted
 * @return std::vector of values
 **/
template <typename T>
inline std::vector<double> value_of_rec(const std::vector<T>& x) {
  size_t size = x.size();
  std::vector<double> result(size);
  for (size_t i = 0; i < size; i++)
    result[i] = value_of_rec(x[i]);
  return result;
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of_rec(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified std::vector.
 * @return Specified std::vector.
 */
template <>
inline std::vector<double> value_of_rec(const std::vector<double>& x) {
  return x;
}

/**
 * Convert a matrix of type T to a matrix of doubles.
 *
 * T must implement value_of_rec. See
 * test/unit/math/fwd/mat/fun/value_of_test.cpp for fvar and var usage.
 *
 * @tparam T Scalar type in matrix
 * @tparam R Rows of matrix
 * @tparam C Columns of matrix
 * @param[in] M Matrix to be converted
 * @return Matrix of values
 **/
template <typename T, int R, int C>
inline Eigen::Matrix<double, R, C> value_of_rec(
    const Eigen::Matrix<T, R, C>& M) {
  Eigen::Matrix<double, R, C> Md(M.rows(), M.cols());
  for (int j = 0; j < M.cols(); j++)
    for (int i = 0; i < M.rows(); i++)
      Md(i, j) = value_of_rec(M(i, j));
  return Md;
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of_rec(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified matrix.
 * @return Specified matrix.
 */
template <int R, int C>
inline typename Eigen::Matrix<double, R, C> value_of_rec(
    const Eigen::Matrix<double, R, C>& x) {
  return x;
}
}  // namespace math
}  // namespace stan

#endif
