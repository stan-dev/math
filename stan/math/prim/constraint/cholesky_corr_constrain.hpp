#ifndef STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <stan/math/prim/fun/tanh.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cstdint>
#include <utility>

namespace stan {
namespace math {

namespace internal {
/**
 * Order of the partial correlations in the unconstrained vector: `row` for
 * `cholesky_corr_constrain` (row 1, then row 2, ...), `column` for
 * `corr_matrix_constrain` and `cov_matrix_constrain_lkj` (column 0, then
 * column 1, ...).
 */
enum class CorrOrder : uint8_t { row = 0, column = 1 };

/**
 * Return the position in the unconstrained vector of the partial
 * correlation for element `(r, c)`, `c < r`, of the Cholesky factor.
 *
 * @tparam Order order of the partial correlations in the vector
 * @param r row
 * @param c column
 * @param K number of rows
 * @return position in the unconstrained vector
 */
template <CorrOrder Order>
inline Eigen::Index cholesky_corr_position(Eigen::Index r, Eigen::Index c,
                                           Eigen::Index K) {
  if constexpr (Order == CorrOrder::column) {
    return c * (K - 1) - (c * (c - 1)) / 2 + (r - c - 1);
  } else {
    return (r * (r - 1)) / 2 + c;
  }
}

/**
 * Return a lazy expression of the coefficients `c` for which the log
 * absolute Jacobian determinant of `cholesky_corr_constrain` is
 * `-sum_n c_n log(cosh(y_n))`, in row order: the coefficient of element
 * `(r, c)` is `r + 1 - c`, with `r` recovered from the position by
 * inverting `cholesky_corr_position<CorrOrder::row>`. The inverse is exact,
 * since at the first position of row `r` the radicand `1 + 8 i` is the
 * perfect square `(2 r - 1)^2`, for which `sqrt` is exact, so `floor` never
 * lands one row too low.
 *
 * @param K number of rows
 * @return coefficients in the order of the unconstrained vector
 */
inline auto cholesky_corr_log_jacobian_coefficients(Eigen::Index K) {
  const Eigen::Index n = (K * (K - 1)) / 2;
  return make_holder(
      [](auto&& i) {
        const auto r = (0.5 * (1.0 + (1.0 + 8.0 * i).sqrt())).floor();
        return (r + 1.0 - (i - 0.5 * r * (r - 1.0))).matrix();
      },
      Eigen::ArrayXd::LinSpaced(n, 0.0, n - 1.0));
}

/**
 * Return a lazy expression of the coefficients `c` for which the log
 * absolute Jacobian determinant of `corr_matrix_constrain` is
 * `-sum_n c_n log(cosh(y_n))`, in column order: the coefficient of every
 * element of column `c` is `K - c` (equation 11 of Lewandowski, Kurowicka
 * and Joe, 2009), with `c` recovered from the position by inverting
 * `cholesky_corr_position<CorrOrder::column>`. The inverse is exact, since
 * at the first position of column `c` the radicand `b^2 - 8 i` is the
 * perfect square `(2 (K - c) - 1)^2`, for which `sqrt` is exact, so `floor`
 * never lands one column too low.
 *
 * @param K number of rows
 * @return coefficients in the order of the unconstrained vector
 */
inline auto corr_matrix_log_jacobian_coefficients(Eigen::Index K) {
  const Eigen::Index n = (K * (K - 1)) / 2;
  const double b = 2.0 * K - 1.0;
  return make_holder(
      [K = static_cast<double>(K), b](auto&& i) {
        return (K - (0.5 * (b - (b * b - 8.0 * i).sqrt())).floor()).matrix();
      },
      Eigen::ArrayXd::LinSpaced(n, 0.0, n - 1.0));
}

/**
 * Return the Cholesky factor of a correlation matrix from unconstrained
 * values `y`, whose partial correlations are `tanh(y)`. Each element is
 * `L(r, c) = tanh(y_c) exp(-a)` with `a` the sum of `log(cosh(y))` over the
 * earlier columns of the row, which has no cancellation.
 *
 * @tparam Order order of the partial correlations in `y`
 * @tparam EigVec type of the vector
 * @tparam LogCoshVec type of the vector of `log(cosh(y))`
 * @param y unconstrained values, `K (K - 1) / 2` of them
 * @param log_cosh_y `log_cosh(y)`, shared with the log Jacobian
 * @param K number of rows
 * @return Cholesky factor of a correlation matrix
 */
template <CorrOrder Order, typename EigVec, typename LogCoshVec,
          require_all_eigen_vector_t<EigVec, LogCoshVec>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_factor(EigVec&& y, LogCoshVec&& log_cosh_y, Eigen::Index K) {
  using T_scalar = value_type_t<EigVec>;
  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> L
      = Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(K, K);
  if (K == 0) {
    return L;
  }
  decltype(auto) y_ref = to_ref(std::forward<EigVec>(y));
  decltype(auto) log_cosh_y_ref = to_ref(std::forward<LogCoshVec>(log_cosh_y));
  L.coeffRef(0, 0) = 1;
  for (Eigen::Index r = 1; r < K; ++r) {
    T_scalar a = 0;
    for (Eigen::Index c = 0; c < r; ++c) {
      const Eigen::Index k = cholesky_corr_position<Order>(r, c, K);
      L.coeffRef(r, c) = tanh(y_ref.coeff(k)) * exp(-a);
      a += log_cosh_y_ref.coeff(k);
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
cholesky_corr_constrain(EigVec&& y, int K) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "constrain size", y.size(),
                   "k_choose_2", k_choose_2);
  decltype(auto) y_ref = to_ref(std::forward<EigVec>(y));
  return internal::cholesky_corr_factor<internal::CorrOrder::row>(
      y_ref, log_cosh(y_ref), K);
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
          require_eigen_col_vector_t<EigVec>* = nullptr,
          require_convertible_t<return_type_t<EigVec>, Lp>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_corr_constrain(EigVec&& y, int K, Lp& lp) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  decltype(auto) y_ref = to_ref(std::forward<EigVec>(y));
  decltype(auto) log_cosh_y = to_ref(log_cosh(y_ref));
  lp -= dot_product(internal::cholesky_corr_log_jacobian_coefficients(K),
                    log_cosh_y);
  return internal::cholesky_corr_factor<internal::CorrOrder::row>(
      y_ref, log_cosh_y, K);
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
