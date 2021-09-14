#ifndef STAN_MATH_PRIM_FUN_COV_MATRIX_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_COV_MATRIX_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the symmetric, positive-definite matrix of dimensions K
 * by K resulting from transforming the specified finite vector of
 * size K plus (K choose 2).
 *
 * <p>See <code>cov_matrix_free()</code> for the inverse transform.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x The vector to convert to a covariance matrix.
 * @param K The number of rows and columns of the resulting
 * covariance matrix.
 * @throws std::invalid_argument if (x.size() != K + (K choose 2)).
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cov_matrix_constrain(const T& x, Eigen::Index K) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;

  Matrix<value_type_t<T>, Dynamic, Dynamic> L(K, K);
  check_size_match("cov_matrix_constrain", "x.size()", x.size(),
                   "K + (K choose 2)", (K * (K + 1)) / 2);
  const auto& x_ref = to_ref(x);
  int i = 0;
  for (Eigen::Index m = 0; m < K; ++m) {
    L.row(m).head(m) = x_ref.segment(i, m);
    i += m;
    L.coeffRef(m, m) = exp(x_ref.coeff(i++));
    L.row(m).tail(K - m - 1).setZero();
  }
  return multiply_lower_tri_self_transpose(L);
}

/**
 * Return the symmetric, positive-definite matrix of dimensions K
 * by K resulting from transforming the specified finite vector of
 * size K plus (K choose 2).
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x The vector to convert to a covariance matrix.
 * @param K The dimensions of the resulting covariance matrix.
 * @param lp Reference
 * @throws std::domain_error if (x.size() != K + (K choose 2)).
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cov_matrix_constrain(const T& x, Eigen::Index K, return_type_t<T>& lp) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using std::log;
  check_size_match("cov_matrix_constrain", "x.size()", x.size(),
                   "K + (K choose 2)", (K * (K + 1)) / 2);
  Matrix<value_type_t<T>, Dynamic, Dynamic> L(K, K);
  const auto& x_ref = to_ref(x);
  int i = 0;
  for (Eigen::Index m = 0; m < K; ++m) {
    L.row(m).head(m) = x_ref.segment(i, m);
    i += m;
    L.coeffRef(m, m) = exp(x_ref.coeff(i++));
    L.row(m).tail(K - m - 1).setZero();
  }
  // Jacobian for complete transform, including exp() above
  lp += (K * LOG_TWO);  // needless constant; want propto
  for (Eigen::Index k = 0; k < K; ++k) {
    lp += (K - k + 1) * log(L.coeff(k, k));  // only +1 because index from 0
  }
  return multiply_lower_tri_self_transpose(L);
}

/**
 * Return the symmetric, positive-definite matrix of dimensions K by K resulting
 * from transforming the specified finite vector of size K plus (K choose 2). If
 * the `Jacobian` parameter is `true`, the log density accumulator is
 * incremented with the log absolute Jacobian determinant of the transform.  All
 * of the transforms are specified with their Jacobians in the *Stan Reference
 * Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column
 * @param x The vector to convert to a covariance matrix
 * @param K The dimensions of the resulting covariance matrix
 * @param[in, out] lp log density accumulator
 * @throws std::domain_error if (x.size() != K + (K choose 2)).
 */
template <bool Jacobian, typename T, require_not_std_vector_t<T>* = nullptr>
inline auto cov_matrix_constrain(const T& x, Eigen::Index K,
                                 return_type_t<T>& lp) {
  if (Jacobian) {
    return cov_matrix_constrain(x, K, lp);
  } else {
    return cov_matrix_constrain(x, K);
  }
}

/**
 * Return the symmetric, positive-definite matrix of dimensions K by K resulting
 * from transforming the specified finite vector of size K plus (K choose 2). If
 * the `Jacobian` parameter is `true`, the log density accumulator is
 * incremented with the log absolute Jacobian determinant of the transform.  All
 * of the transforms are specified with their Jacobians in the *Stan Reference
 * Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param x The vector to convert to a covariance matrix
 * @param K The dimensions of the resulting covariance matrix
 * @param[in, out] lp log density accumulator
 * @throws std::domain_error if (x.size() != K + (K choose 2)).
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto cov_matrix_constrain(const T& x, Eigen::Index K,
                                 return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(x, [&lp, K](auto&& v) {
    return cov_matrix_constrain<Jacobian>(v, K, lp);
  });
}

}  // namespace math
}  // namespace stan

#endif
