#ifndef STAN_MATH_PRIM_FUN_REP_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_REP_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Implementation of rep_matrix returning an Eigen matrix with scalar
 * type equal to the input scalar type.
 * @tparam Ret An Eigen type.
 * @tparam T A Scalar type.
 * @param x A Scalar whose values are propogated to all values in the return
 * matrix.
 * @param m Number or rows.
 * @param n Number of columns.
 */
template <typename Ret, typename T,
          require_eigen_matrix_dynamic_vt<is_stan_scalar, Ret>* = nullptr,
          require_stan_scalar_t<T>* = nullptr>
inline auto rep_matrix(const T& x, int m, int n) {
  check_nonnegative("rep_matrix", "rows", m);
  check_nonnegative("rep_matrix", "cols", n);
  return Ret::Constant(m, n, x);
}

/**
 * Default Implementation of rep_matrix returning an Eigen matrix with scalar
 * type equal to the input scalar type.
 * @tparam T A Scalar type.
 * @param x A Scalar whose values are propogated to all values in the return
 * matrix.
 * @param m Number or rows.
 * @param n Number of columns.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline auto rep_matrix(const T& x, int m, int n) {
  return rep_matrix<
      Eigen::Matrix<return_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>>(x, m, n);
}

/**
 * Implementation of rep_matrix returning an Eigen matrix from an Eigen
 * vector.
 * @tparam Vec An Eigen vector.
 * @param x An Eigen vector. For Row vectors the values are replacated rowwise.
 * and for column vectors the values are repliacated colwise.
 * @param n Number of rows or columns.
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
inline auto rep_matrix(const Vec& x, int n) {
  if (is_eigen_row_vector<Vec>::value) {
    check_nonnegative("rep_matrix", "rows", n);
    return x.replicate(n, 1);
  } else {
    check_nonnegative("rep_matrix", "cols", n);
    return x.replicate(1, n);
  }
}

}  // namespace math
}  // namespace stan

#endif
