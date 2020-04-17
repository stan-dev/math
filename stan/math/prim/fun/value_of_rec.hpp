#ifndef STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_PRIM_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <complex>
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
 * @tparam T Type of scalar.
 * @param x Scalar to convert to double.
 * @return Value of scalar cast to a double.
 */
template <typename T, require_floating_point_t<T>* = nullptr>
inline auto value_of_rec(T x) {
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
 * @param x Specified value.
 * @return Specified value.
 */
template <typename T, require_integral_t<T>* = nullptr>
inline auto value_of_rec(T x) {
  return x;
}

/**
 * Recursively apply value-of to the parts of the argument.
 *
 * @tparam T value type of argument
 * @param[in] x argument
 * @return real complex value of argument
 */
template <typename ComplexT, require_complex_t<ComplexT>* = nullptr,
          require_not_vt_arithmetic<ComplexT>* = nullptr>
inline auto value_of_rec(ComplexT&& x) {
  return std::complex<double>{value_of_rec(x.real()), value_of_rec(x.imag())};
}

/**
 * Recursively apply value-of to the parts of the argument.
 *
 * @tparam ComplexT value type of argument
 * @param[in] x argument
 * @return real complex value of argument
 */
template <typename ComplexT, require_complex_t<ComplexT>* = nullptr,
          require_vt_arithmetic<ComplexT>* = nullptr>
inline decltype(auto) value_of_rec(ComplexT&& x) {
  return std::forward<ComplexT>(x);
}

/**
 * Convert a std::vector of type T to a std::vector of doubles.
 *
 * T must implement value_of_rec. See
 * test/math/fwd/fun/value_of_rec.cpp for fvar and var usage.
 *
 * @tparam T Scalar type in std::vector
 * @param[in] x std::vector to be converted
 * @return std::vector of values
 **/
template <typename Vec, require_std_vector_t<Vec>* = nullptr,
          require_not_vt_arithmetic<Vec>* = nullptr,
          require_not_vt_var<Vec>* = nullptr>
inline auto value_of_rec(Vec&& x) {
  size_t x_size = x.size();
  std::vector<double> result(x_size);
  for (size_t i = 0; i < x_size; i++) {
    result[i] = value_of_rec(x[i]);
  }
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
template <typename Vec,
          require_std_vector_vt<std::is_arithmetic, Vec>* = nullptr>
inline decltype(auto) value_of_rec(Vec&& x) {
  return std::forward<Vec>(x);
}

/**
 * Convert a matrix of type T to a matrix of doubles.
 *
 * T must implement value_of_rec. See
 * test/unit/math/fwd/fun/value_of_test.cpp for fvar and var usage.
 *
 * @tparam T Type of matrix
 * @param[in] M Matrix to be converted
 * @return Matrix of values
 **/
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_arithmetic<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline auto value_of_rec(EigMat&& M) {
  return M.unaryExpr([](auto&& x) { return value_of_rec(x); });
}

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of_rec(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @tparam T Type of matrix.
 * @param x Specified matrix.
 * @return Specified matrix.
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_vt_arithmetic<EigMat>* = nullptr>
inline decltype(auto) value_of_rec(EigMat&& M) {
  return std::forward<EigMat>(M);
}

}  // namespace math
}  // namespace stan

#endif
