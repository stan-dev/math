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
 * This specialization handles all floating point types and passes the
 * input to the output.
 *
 * @tparam T type that is or contains a floating point value.
 * @param x scalar to convert to double
 * @return value of scalar cast to double
 */
template <typename T, require_vt_floating_point<T>* = nullptr>
inline decltype(auto) value_of(T&& x) {
  return std::forward<T>(x);
}

/**
 * Specialization to convert integral types to double.
 *
 * @param x value
 * @return input value
 */
template <typename T, require_integral_t<T>* = nullptr>
inline auto value_of(T x) {
  return static_cast<double>(x);
}

/**
 * Return the inner value of the specified complex type.
 * @tparam ComplexT The complex type.
 * @param x A complex object.
 * @return A complex with `value_of` applied to the `real` and `imag` parts.
 */
template <typename ComplexT, require_complex_t<ComplexT>* = nullptr,
          require_not_vt_arithmetic<ComplexT>* = nullptr>
inline auto value_of(ComplexT&& x) {
  using complex_ret = std::complex<partials_type_t<value_type_t<ComplexT>>>;
  return complex_ret{value_of(x.real()), value_of(x.imag())};
}

/**
 * Apply `value_of` to all elements of a `std::vector`
 *
 * @tparam Vec A vector containing values with a defined `value_of`.
 * @param[in] x std::vector to be converted
 * @return std::vector with `value_of` applied to each element.
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
 * Specialization to convert vector of integral types to double.
 * @tparam Vec Standard vector holding integral type.
 * @param x Standard vector.
 */
template <typename Vec, require_std_vector_vt<std::is_integral, Vec>* = nullptr>
inline decltype(auto) value_of(Vec&& x) {
  return std::forward<Vec>(x);
}

/**
 * Apply `value_of` to each element of an Eigen type.
 *
 * Scalar type of matrix must impliment `value_of`. See
 * test/math/fwd/fun/value_of.cpp for fvar and var usage.
 *
 * @tparam EigMat type derived from `EigenBase` whose scalar type is neither
 * var or arithmetic.
 *
 * @param[in] M Matrix to apply over.
 * @return Matrix of values
 **/
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr,
          require_not_vt_arithmetic<EigMat>* = nullptr>
inline auto value_of(EigMat&& M) {
  using eig_mat = std::decay_t<EigMat>;
  using eig_partial = partials_type_t<value_type_t<EigMat>>;
  constexpr Eigen::Index R = eig_mat::RowsAtCompileTime;
  constexpr Eigen::Index C = eig_mat::ColsAtCompileTime;
  Eigen::Matrix<eig_partial, R, C> Md = M.unaryExpr([&](auto&& x) {
    return value_of(x);
  });
  return Md;
}

/**
 * Specialization to convert integral matrices types to double matrices.
 *
 * @param x Specified matrix.
 * @return Specified matrix.
 */
template <typename EigMat,
          require_eigen_vt<std::is_integral, EigMat>* = nullptr>
inline auto value_of(EigMat&& x) {
  return x.template cast<double>();
}

}  // namespace math
}  // namespace stan

#endif
