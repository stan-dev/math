#ifndef STAN_MATH_REV_FUN_REP_MATRIX_HPP
#define STAN_MATH_REV_FUN_REP_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/rep_matrix.hpp>

namespace stan {
namespace math {

/**
 * Impl of rep_matrix returning an `var_value<Eigen::Matrix>` with a double
 * scalar type.
 * @tparam Ret A `var_value` with inner Eigen type.
 * @tparam T A Scalar type.
 * @param x A Scalar whose values are propogated to all values in the return
 * matrix.
 * @param m Number or rows.
 * @param n Number of columns.
 */
template <typename Ret, typename T, require_var_matrix_t<Ret>* = nullptr,
          require_arithmetic_t<T>* = nullptr>
inline auto rep_matrix(const T& x, int m, int n) {
  check_nonnegative("rep_matrix", "rows", m);
  check_nonnegative("rep_matrix", "cols", n);
  return Ret(value_type_t<Ret>::Constant(m, n, x));
}

/**
 * Impl of rep_matrix returning an `var_value<Eigen::Matrix>` with a var scalar
 * type.
 * @tparam Ret A `var_value` with inner Eigen type.
 * @tparam T A Scalar type.
 * @param x A Scalar whose values are propogated to all values in the return
 * matrix.
 * @param m Number or rows.
 * @param n Number of columns.
 */
template <typename Ret, typename T, require_var_matrix_t<Ret>* = nullptr,
          require_var_t<T>* = nullptr>
inline auto rep_matrix(const T& x, int m, int n) {
  check_nonnegative("rep_matrix", "rows", m);
  check_nonnegative("rep_matrix", "cols", n);
  return make_callback_var(
      value_type_t<Ret>::Constant(m, n, x.val()),
      [x](auto& rep) mutable { x.adj() += rep.adj().sum(); });
}

/**
 * Impl of rep_matrix returning an `var_value<Eigen::Matrix>` from an Eigen
 * matrix with Arithmetic scalar.
 * @tparam Ret A `var_value` with inner Eigen dynamic matrix type.
 * @tparam Vec An Eigen vector with Arithmetic scalar type.
 * @param x An Eigen vector. For Row vectors the values are replacated rowwise
 * and for column vectors the values are repliacated colwise.
 * @param n Number of rows or columns.
 */
template <typename Ret, typename Vec, require_var_matrix_t<Ret>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, Vec>* = nullptr>
inline auto rep_matrix(const Vec& x, int n) {
  if (is_eigen_row_vector<Vec>::value) {
    check_nonnegative("rep_matrix", "rows", n);
    return Ret(x.replicate(n, 1));
  } else {
    check_nonnegative("rep_matrix", "cols", n);
    return Ret(x.replicate(1, n));
  }
}

/**
 * Impl of rep_matrix returning a `var_value<Eigen::Matrix>` from a `var_value`
 * with an inner Eigen vector type.
 * @tparam Ret A `var_value` with inner Eigen dynamic matrix type.
 * @tparam Vec A `var_value` with an inner Eigen vector type.
 * @param x A `var_value` with inner Eigen vector type. For Row vectors the
 * values are replacated rowwise and for column vectors the values are
 * repliacated colwise.
 * @param n Number of rows or columns.
 */
template <typename Ret, typename Vec, require_var_matrix_t<Ret>* = nullptr,
          require_vector_st<is_var, Vec>* = nullptr>
inline auto rep_matrix(const Vec& x, int n) {
  if (is_row_vector<Vec>::value) {
    check_nonnegative("rep_matrix", "rows", n);
    return make_callback_var(x.val().replicate(n, 1), [x](auto& rep) mutable {
      x.adj() += rep.adj().rowwise().sum();
    });
  } else {
    check_nonnegative("rep_matrix", "cols", n);
    return make_callback_var(x.val().replicate(1, n), [x](auto& rep) mutable {
      x.adj() += rep.adj().colwise().sum();
    });
  }
}

}  // namespace math
}  // namespace stan

#endif
