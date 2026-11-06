#ifndef STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <stan/math/prim/fun/tanh.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <utility>

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

/**
 * Return `tanh(y)` and `sech(y) = 1 / cosh(y)` from one call of `expm1` or
 * `exp`, with a relative error of a few ulp each.
 *
 * For `|y| < 1`, with `u = expm1(|y|)`, `cosh = 1 + u^2 / (2 (1 + u))` and
 * `sinh = u (u + 2) / (2 (1 + u))`. For `|y| >= 1`, with `e = exp(-|y|)`,
 * `tanh = (1 - e^2) / (1 + e^2)` and `sech = 2 e / (1 + e^2)`.
 *
 * @param y argument
 * @return pair of `tanh(y)` and `sech(y)`
 */
inline std::pair<double, double> tanh_and_sech(double y) {
  const double a = std::fabs(y);
  double tanh_a;
  double sech;
  if (a < 1.0) {
    const double u = std::expm1(a);
    const double cosh = 1.0 + u * u / (2.0 * (1.0 + u));
    tanh_a = u * (u + 2.0) / (2.0 * (1.0 + u)) / cosh;
    sech = 1.0 / cosh;
  } else {
    const double e = std::exp(-a);
    const double e2 = e * e;
    tanh_a = (1.0 - e2) / (1.0 + e2);
    sech = 2.0 * e / (1.0 + e2);
  }
  return {std::copysign(tanh_a, y), sech};
}

/**
 * Return the position in the unconstrained vector of the partial
 * correlation for element `(r, c)`, `c < r`, of the Cholesky factor.
 *
 * @tparam ColumnOrder `false` for the row order of
 * `cholesky_corr_constrain` (row 1, then row 2, ...), `true` for the column
 * order of `corr_matrix_constrain` (column 0, then column 1, ...)
 * @param r row
 * @param c column
 * @param K number of rows
 * @return position in the unconstrained vector
 */
template <bool ColumnOrder>
inline Eigen::Index cholesky_corr_position(Eigen::Index r, Eigen::Index c,
                                           Eigen::Index K) {
  if constexpr (ColumnOrder) {
    return c * (K - 1) - (c * (c - 1)) / 2 + (r - c - 1);
  } else {
    return (r * (r - 1)) / 2 + c;
  }
}

/**
 * Return the Cholesky factor of a correlation matrix from unconstrained
 * values `y`, whose partial correlations are `tanh(y)`.
 *
 * In row `r`, `1 - sum` of the squares of the elements before column `c`
 * is `prod_{m < c} (1 - tanh(y_m)^2) = exp(-2 a)` with
 * `a = sum_{m < c} log(cosh(y_m))`, so
 * `L(r, c) = tanh(y_c) exp(-a)` and `L(r, r) = exp(-a)` over the whole
 * row. Unlike `sqrt(1 - sum)`, this has no subtraction: the elements keep
 * a small relative error when partial correlations are close to ±1 or a
 * row is long, and the diagonal stays positive unless `a` exceeds about
 * 745.
 *
 * @tparam ColumnOrder order of the partial correlations in `y`, see
 * `cholesky_corr_position`
 * @tparam EigVec type of the vector
 * @param y unconstrained values, `K (K - 1) / 2` of them
 * @param K number of rows
 * @return Cholesky factor of a correlation matrix
 */
template <bool ColumnOrder, typename EigVec,
          require_eigen_vector_t<EigVec>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_factor(const EigVec& y, Eigen::Index K) {
  using T_scalar = value_type_t<EigVec>;
  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> L
      = Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(K, K);
  if (K == 0) {
    return L;
  }
  const auto& y_ref = to_ref(y);
  L.coeffRef(0, 0) = 1;
  for (Eigen::Index r = 1; r < K; ++r) {
    T_scalar a = 0;
    for (Eigen::Index c = 0; c < r; ++c) {
      const auto& y_rc
          = y_ref.coeff(cholesky_corr_position<ColumnOrder>(r, c, K));
      L.coeffRef(r, c) = tanh(y_rc) * exp(-a);
      a += log_cosh(y_rc);
    }
    L.coeffRef(r, r) = exp(-a);
  }
  return L;
}
}  // namespace internal

/**
 * Return the Cholesky factor of a `K x K` correlation matrix from the
 * `K (K - 1) / 2` unconstrained values `y`, whose partial correlations are
 * `tanh(y)` in row order.
 *
 * @tparam EigVec type of the vector
 * @param y unconstrained values
 * @param K number of rows
 * @return Cholesky factor of a correlation matrix
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_constrain(const EigVec& y, int K) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "constrain size", y.size(),
                   "k_choose_2", k_choose_2);
  return internal::cholesky_corr_factor<false>(y, K);
}

/**
 * Return the Cholesky factor of a `K x K` correlation matrix from the
 * `K (K - 1) / 2` unconstrained values `y`, and add the log absolute
 * Jacobian determinant of the transform to `lp`.
 *
 * @tparam EigVec type of the vector
 * @tparam Lp scalar type of the log density accumulator
 * @param y unconstrained values
 * @param K number of rows
 * @param[in,out] lp log density accumulator
 * @return Cholesky factor of a correlation matrix
 */
template <typename EigVec, typename Lp,
          require_eigen_vector_t<EigVec>* = nullptr,
          require_convertible_t<return_type_t<EigVec>, Lp>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_constrain(const EigVec& y, int K, Lp& lp) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  if (K == 0) {
    return Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>(
        0, 0);
  }
  const auto& y_ref = to_ref(y);
  lp -= dot_product(internal::cholesky_corr_log_jacobian_coefficients(K),
                    log_cosh(y_ref));
  return internal::cholesky_corr_factor<false>(y_ref, K);
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
