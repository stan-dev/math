#ifndef STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/constraint/corr_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
/**
 * Return the coefficients `c` for which the log absolute Jacobian
 * determinant of `cholesky_corr_constrain(y, K)` is
 * `-sum_n c_n log(cosh(y_n))`.
 *
 * With `z = tanh(y)`, row `i` holds `z_m sqrt(prod_{l < m} (1 - z_l^2))` in
 * column `m < i`. The Jacobian has a factor `1 - z_m^2` for each element
 * and a factor `sqrt(prod_{l < j} (1 - z_l^2))` for each column
 * `0 < j < i`. With `log(1 - z_m^2) = -2 log(cosh(y_m))`, the coefficient
 * of element `m` of row `i` is `2 + (i - 1 - m) = i + 1 - m`.
 *
 * @param K number of rows of the Cholesky factor
 * @return coefficients in the order of the unconstrained vector
 */
inline Eigen::VectorXd cholesky_corr_log_jacobian_coefficients(int K) {
  Eigen::VectorXd c((K * (K - 1)) / 2);
  int k = 0;
  for (int i = 1; i < K; ++i) {
    for (int m = 0; m < i; ++m) {
      c.coeffRef(k++) = i + 1 - m;
    }
  }
  return c;
}
}  // namespace internal

template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_constrain(const EigVec& y, int K) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::sqrt;
  using T_scalar = value_type_t<EigVec>;
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "constrain size", y.size(),
                   "k_choose_2", k_choose_2);
  Matrix<T_scalar, Dynamic, 1> z = corr_constrain(y);
  Matrix<T_scalar, Dynamic, Dynamic> x(K, K);
  if (K == 0) {
    return x;
  }
  x.setZero();
  x.coeffRef(0, 0) = 1;
  int k = 0;
  for (int i = 1; i < K; ++i) {
    x.coeffRef(i, 0) = z.coeff(k++);
    T_scalar sum_sqs = square(x.coeff(i, 0));
    for (int j = 1; j < i; ++j) {
      x.coeffRef(i, j) = z.coeff(k++) * sqrt(1.0 - sum_sqs);
      sum_sqs += square(x.coeff(i, j));
    }
    x.coeffRef(i, i) = sqrt(1.0 - sum_sqs);
  }
  return x;
}

// FIXME to match above after debugged
template <typename EigVec, typename Lp,
          require_eigen_vector_t<EigVec>* = nullptr,
          require_convertible_t<return_type_t<EigVec>, Lp>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_constrain(const EigVec& y, int K, Lp& lp) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::sqrt;
  using T_scalar = value_type_t<EigVec>;
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  Matrix<T_scalar, Dynamic, 1> z = corr_constrain(y);
  Matrix<T_scalar, Dynamic, Dynamic> x(K, K);
  if (K == 0) {
    return x;
  }
  lp -= dot_product(internal::cholesky_corr_log_jacobian_coefficients(K),
                    log_cosh(y));
  x.setZero();
  x.coeffRef(0, 0) = 1;
  int k = 0;
  for (int i = 1; i < K; ++i) {
    x.coeffRef(i, 0) = z.coeff(k++);
    T_scalar sum_sqs = square(x.coeff(i, 0));
    for (int j = 1; j < i; ++j) {
      x.coeffRef(i, j) = z.coeff(k++) * sqrt(1.0 - sum_sqs);
      sum_sqs += square(x.coeff(i, j));
    }
    x.coeffRef(i, i) = sqrt(1.0 - sum_sqs);
  }
  return x;
}

/**
 * Return The cholesky of a `KxK` correlation matrix.
 * This overload handles looping over the elements of a standard vector.
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param y Linearly Serialized vector of size `(K * (K - 1))/2` holding the
 *  column major order elements of the lower triangurlar
 * @param K The size of the matrix to return
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto cholesky_corr_constrain(T&& y, int K) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [K](auto&& v) {
    return cholesky_corr_constrain(std::forward<decltype(v)>(v), K);
  });
}

/**
 * Return The cholesky of a `KxK` correlation matrix.
 * This overload handles looping over the elements of a standard vector.
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @tparam Lp Scalar type for the lp argument. The scalar type of T should be
 * convertible to this.
 * @param y Linearly Serialized vector of size `(K * (K - 1))/2` holding the
 *  column major order elements of the lower triangurlar
 * @param K The size of the matrix to return
 * @param[in,out] lp log density accumulator
 */
template <typename T, typename Lp, require_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto cholesky_corr_constrain(T&& y, int K, Lp& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(y), [&lp, K](auto&& v) {
    return cholesky_corr_constrain(std::forward<decltype(v)>(v), K, lp);
  });
}

/**
 * Return The cholesky of a `KxK` correlation matrix. If the `Jacobian`
 * parameter is `true`, the log density accumulator is incremented with the log
 * absolute Jacobian determinant of the transform.  All of the transforms are
 * specified with their Jacobians in the *Stan Reference Manual* chapter
 * Constraint Transforms.
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column, or a standard vector thereof
 * @tparam Lp A scalar type for the lp argument. The scalar type of T should be
 * convertible to this.
 * @param y Linearly Serialized vector of size `(K * (K - 1))/2` holding the
 *  column major order elements of the lower triangurlar
 * @param K The size of the matrix to return
 * @param[in,out] lp log density accumulator
 */
template <bool Jacobian, typename T, typename Lp,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto cholesky_corr_constrain(T&& y, int K, Lp& lp) {
  if constexpr (Jacobian) {
    return cholesky_corr_constrain(std::forward<T>(y), K, lp);
  } else {
    return cholesky_corr_constrain(std::forward<T>(y), K);
  }
}

}  // namespace math
}  // namespace stan
#endif
