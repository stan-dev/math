#ifndef STAN_MATH_REV_FUNCTOR_APPLY_SCALAR_BINARY_HPP
#define STAN_MATH_REV_FUNCTOR_APPLY_SCALAR_BINARY_HPP

#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/num_elements.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Specialisation for use with combinations of
 * `Eigen::Matrix` and `var_value<Eigen::Matrix>` inputs.
 * Eigen's binaryExpr framework is used for more efficient indexing of both row-
 * and column-major inputs  without separate loops.
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First Matrix input to which operation is applied.
 * @param y Second Matrix input to which operation is applied.
 * @param f functor to apply to Matrix inputs.
 * @return `var_value<Matrix>` with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_any_var_matrix_t<T1, T2>* = nullptr,
          require_all_matrix_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  check_matching_dims("Binary function", "x", x, "y", y);
  return f(x, y);
}

/**
 * Specialisation for use with one `var_value<Eigen vector>` (row or column) and
 * a one-dimensional std::vector of integer types
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Matrix input to which operation is applied.
 * @param y Integer std::vector input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return var_value<Eigen> object with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_var_matrix_t<T1>* = nullptr,
          require_std_vector_vt<std::is_integral, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  return f(x, y);
}

/**
 * Specialisation for use with a one-dimensional std::vector of integer types
 * and one `var_value<Eigen vector>` (row or column).
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Integer std::vector input to which operation is applied.
 * @param y Eigen input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_std_vector_vt<std::is_integral, T1>* = nullptr,
          require_var_matrix_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  return f(x, y);
}

/**
 * Specialisation for use with one `var_value<Matrix>` and
 * a two-dimensional std::vector of integer types
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x var with Eigen matrix inner type to which operation is applied.
 * @param y Nested integer std::vector input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_var_matrix_t<T1>* = nullptr,
          require_std_vector_vt<is_std_vector, T2>* = nullptr,
          require_std_vector_st<std::is_integral, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return f(x, y);
}

/**
 * Specialisation for use with a two-dimensional std::vector of integer types
 * and one `var_value<Matrix>`.
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Nested integer std::vector input to which operation is applied.
 * @param y var value with inner Eigen matrix input to which operation is
 * applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_std_vector_vt<is_std_vector, T1>* = nullptr,
          require_std_vector_st<std::is_integral, T1>* = nullptr,
          require_var_matrix_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return f(x, y);
}

/**
 * Specialisation for use when the first input is an `var_value<Eigen> type and
 * the second is a scalar.
 *
 * @tparam T1 Type of `var_value<Matrix>` object to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Matrix input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to var matrix and scalar inputs.
 * @return `var_value<Matrix> object with result of applying functor to inputs.
 *
 * Note: The return expresssion needs to be evaluated, otherwise the captured
 *         function and scalar fall out of scope.
 */
template <typename T1, typename T2, typename F,
          require_var_matrix_t<T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return f(x, y);
}

/**
 * Specialisation for use when the first input is an scalar and the second is
 * an `var_value<Eigen>`.
 *
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 var value with inner Eigen type to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y var matrix input to which operation is applied.
 * @param f functor to apply to Eigen and scalar inputs.
 * @return var value with inner Eigen type with result of applying functor to
 * inputs.
 *
 * Note: The return expresssion needs to be evaluated, otherwise the captured
 *         function and scalar fall out of scope.
 */
template <typename T1, typename T2, typename F,
          require_stan_scalar_t<T1>* = nullptr,
          require_var_matrix_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return f(x, y);
}

/**
 * Specialisation for use when the first input is a nested std::vector and the
 * second is a scalar. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam T1 Type of std::vector to which functor is applied.
 * @tparam T2 Type of scalar to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x std::vector input to which operation is applied.
 * @param y Scalar input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_std_vector_t<T1>* = nullptr,
          require_var_matrix_t<value_type_t<T1>>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_binary(x[0], y, f))>;
  size_t x_size = x.size();
  std::vector<T_return> result(x_size);
  for (size_t i = 0; i < x_size; ++i) {
    result[i] = apply_scalar_binary(x[i], y, f);
  }
  return result;
}

/**
 * Specialisation for use when the first input is a scalar and the second is a
 * nested std::vector. The returned scalar type is deduced to allow for cases
 * where the input and return scalar types differ (e.g., functions implicitly
 * promoting integers).
 *
 * @tparam T1 Type of scalar to which functor is applied.
 * @tparam T2 Type of std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Scalar input to which operation is applied.
 * @param y std::vector input to which operation is applied.
 * @param f functor to apply to inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_stan_scalar_t<T1>* = nullptr,
          require_std_vector_t<T2>* = nullptr,
          require_var_matrix_t<value_type_t<T2>>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  using T_return = plain_type_t<decltype(apply_scalar_binary(x, y[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(x, y[i], f);
  }
  return result;
}

/**
 * Specialisation for use with two nested containers (std::vectors).
 * The returned scalar type is deduced to allow for cases where the input and
 * return scalar types differ (e.g., functions implicitly promoting
 * integers).
 *
 * @tparam T1 Type of first std::vector to which functor is applied.
 * @tparam T2 Type of second std::vector to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x First std::vector input to which operation is applied.
 * @param y Second std::vector input to which operation is applied.
 * @param f functor to apply to std::vector inputs.
 * @return std::vector with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_any_var_matrix_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(const std::vector<T1>& x,
                                const std::vector<T2>& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  using T_return = plain_type_t<decltype(apply_scalar_binary(x[0], y[0], f))>;
  size_t y_size = y.size();
  std::vector<T_return> result(y_size);
  for (size_t i = 0; i < y_size; ++i) {
    result[i] = apply_scalar_binary(x[i], y[i], f);
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
