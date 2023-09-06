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
 * Specialization for use with combinations of
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
 * Specialization for use with one `var_value<Eigen vector>` (row or column) and
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
          require_any_var_matrix_t<T1, T2>* = nullptr,
          require_any_std_vector_vt<std::is_integral, T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  check_matching_sizes("Binary function", "x", x, "y", y);
  return f(x, y);
}

/**
 * Specialization for use with a two-dimensional std::vector of integer types
 * and one `var_value<Matrix>`.
 *
 * @tparam T1 Type of first argument to which functor is applied.
 * @tparam T2 Type of second argument to which functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Either a var matrix or nested integer std::vector input to which
 * operation is applied.
 * @param x Either a var matrix or nested integer std::vector input to which
 * operation is applied.
 * @param f functor to apply to inputs.
 * @return Eigen object with result of applying functor to inputs.
 */
template <typename T1, typename T2, typename F,
          require_any_std_vector_vt<is_std_vector, T1, T2>* = nullptr,
          require_any_std_vector_st<std::is_integral, T1, T2>* = nullptr,
          require_any_var_matrix_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return f(x, y);
}

/**
 * Specialization for use when the one input is an `var_value<Eigen> type and
 * the other is a scalar.
 *
 * @tparam T1 Type of either `var_value<Matrix>` or scalar object to which
 * functor is applied.
 * @tparam T2 Type of either `var_value<Matrix>` or scalar object to which
 * functor is applied.
 * @tparam F Type of functor to apply.
 * @param x Matrix or Scalar input to which operation is applied.
 * @param x Matrix or Scalar input to which operation is applied.
 * @param f functor to apply to var matrix and scalar inputs.
 * @return `var_value<Matrix> object with result of applying functor to inputs.
 *
 */
template <typename T1, typename T2, typename F,
          require_any_stan_scalar_t<T1, T2>* = nullptr,
          require_any_var_matrix_t<T1, T2>* = nullptr>
inline auto apply_scalar_binary(const T1& x, const T2& y, const F& f) {
  return f(x, y);
}

}  // namespace math
}  // namespace stan
#endif
