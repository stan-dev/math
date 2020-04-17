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
template <typename T, require_t<std::is_floating_point<T>>* = nullptr, require_not_same_t<double, T>* = nullptr>
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
template <typename T, require_same_t<double, T>* = nullptr>
inline decltype(auto) value_of(T&& x) {
  return std::forward<T>(x);
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
template <typename T, require_same_t<int, T>* = nullptr>
inline decltype(auto) value_of(T&& x) { return x; }

template <typename T, require_complex_t<T>* = nullptr, require_not_vt_arithmetic<T>* = nullptr>
inline std::complex<partials_type_t<T>> value_of(T&& x) {
  return {value_of(x.real()), value_of(x.imag())};
}

template <typename T, require_complex_t<T>* = nullptr, require_vt_arithmetic<T>* = nullptr>
inline decltype(auto) value_of(T&& x) {
  return std::forward<T>(x); 
}

/**
 * Convert a std::vector of type T to a std::vector of
 * child_type<T>::type.
 *
 * @tparam T Scalar type in std::vector
 * @param[in] x std::vector to be converted
 * @return std::vector of values
 **/
template <typename Vec, require_std_vector_t<Vec>* = nullptr,
 require_not_vt_arithmetic<Vec>* = nullptr>
inline auto value_of(Vec&& x) {
  const size_t x_size = x.size();
  std::vector<partials_type_t<value_type_t<Vec>>> result(x_size);
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
template <typename Vec, require_std_vector_vt<std::is_arithmetic, Vec>* = nullptr>
inline decltype(auto) value_of(Vec&& x) {
  return std::forward<Vec>(x);
}
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
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
  require_not_vt_var<EigMat>* = nullptr,
  require_not_vt_arithmetic<EigMat>* = nullptr>
inline auto value_of(EigMat&& M) {
  using ref_inner = const typename std::decay_t<EigMat>::PlainObject;
  Eigen::Matrix<partials_type_t<value_type_t<EigMat>>,
   std::decay_t<EigMat>::RowsAtCompileTime, std::decay_t<EigMat>::ColsAtCompileTime> Md(M.rows(), M.cols());
 const Eigen::Ref<ref_inner, Eigen::Aligned16, Eigen::Stride<0,0>>& mat = M;
  for (int j = 0; j < mat.cols(); j++) {
    for (int i = 0; i < mat.rows(); i++) {
      Md(i, j) = value_of(mat.coeffRef(i, j));
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
template <typename EigMat, require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline decltype(auto) value_of(EigMat&& x) {
  return std::forward<EigMat>(x);
}


}  // namespace math
}  // namespace stan

#endif
