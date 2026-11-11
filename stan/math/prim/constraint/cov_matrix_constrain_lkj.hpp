#ifndef STAN_MATH_PRIM_CONSTRAINT_COV_MATRIX_CONSTRAIN_LKJ_HPP
#define STAN_MATH_PRIM_CONSTRAINT_COV_MATRIX_CONSTRAIN_LKJ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/cholesky_corr_constrain.hpp>
#include <stan/math/prim/constraint/positive_constrain.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/diag_pre_multiply.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the `k x k` covariance matrix `(diag(sds) L) (diag(sds) L)^T`, with
 * `L` the Cholesky factor of the correlation matrix of the first
 * `k choose 2` unconstrained values and `sds = exp(x.tail(k))`. See
 * `corr_matrix_constrain` for the correlation transform.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time column)
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations.
 * @param k Dimensionality of returned covariance matrix.
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cov_matrix_constrain_lkj(T&& x, size_t k) {
  size_t k_choose_2 = (k * (k - 1)) / 2;
  decltype(auto) x_ref = to_ref(std::forward<T>(x));
  check_size_match("cov_matrix_constrain_lkj", "x.size()", x_ref.size(),
                   "k_choose_2 + k", k_choose_2 + k);
  auto&& cpcs = x_ref.head(k_choose_2);
  return multiply_lower_tri_self_transpose(diag_pre_multiply(
      to_ref(positive_constrain(x_ref.tail(k))),
      internal::cholesky_corr_factor<internal::CorrOrder::column>(
          cpcs, log_cosh(cpcs), k)));
}

/**
 * Return the covariance matrix as `cov_matrix_constrain_lkj(x, k)` and add
 * the log absolute Jacobian determinant to `lp`. The Jacobian is the one of
 * `corr_matrix_constrain` plus `k log 2 + (k + 1) sum(x.tail(k))`.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time column)
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations.
 * @param k Dimensionality of returned covariance matrix.
 * @param lp Log probability reference to increment.
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cov_matrix_constrain_lkj(T&& x, size_t k, return_type_t<T>& lp) {
  size_t k_choose_2 = (k * (k - 1)) / 2;
  decltype(auto) x_ref = to_ref(std::forward<T>(x));
  check_size_match("cov_matrix_constrain_lkj", "x.size()", x_ref.size(),
                   "k_choose_2 + k", k_choose_2 + k);
  auto&& cpcs = x_ref.head(k_choose_2);
  auto&& log_sds = x_ref.tail(k);
  decltype(auto) log_cosh_cpcs = to_ref(log_cosh(cpcs));
  lp -= dot_product(internal::corr_matrix_log_jacobian_coefficients(k),
                    log_cosh_cpcs);
  lp += (k + 1.0) * sum(log_sds) + k * LOG_TWO;
  return multiply_lower_tri_self_transpose(diag_pre_multiply(
      to_ref(positive_constrain(log_sds)),
      internal::cholesky_corr_factor<internal::CorrOrder::column>(
          cpcs, log_cosh_cpcs, k)));
}

/**
 * Return the covariance matrix of the specified dimensionality derived from
 * constraining the specified vector of unconstrained values. If the `Jacobian`
 * parameter is `true`, the log density accumulator is incremented with the log
 * absolute Jacobian determinant of the transform.  All of the transforms are
 * specified with their Jacobians in the *Stan Reference Manual* chapter
 * Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time rows or
 * columns equal to 1
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations
 * @param k Dimensionality of returned covariance matrix
 * @param[in, out] lp log density accumulator
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <bool Jacobian, typename T, require_not_std_vector_t<T>* = nullptr>
inline auto cov_matrix_constrain_lkj(T&& x, size_t k, return_type_t<T>& lp) {
  if constexpr (Jacobian) {
    return cov_matrix_constrain_lkj(std::forward<T>(x), k, lp);
  } else {
    return cov_matrix_constrain_lkj(std::forward<T>(x), k);
  }
}

/**
 * Return the covariance matrix of the specified dimensionality derived from
 * constraining the specified vector of unconstrained values. If the `Jacobian`
 * parameter is `true`, the log density accumulator is incremented with the log
 * absolute Jacobian determinant of the transform.  All of the transforms are
 * specified with their Jacobians in the *Stan Reference Manual* chapter
 * Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time rows or columns equal to 1
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations
 * @param k Dimensionality of returned covariance matrix
 * @param[in, out] lp log density accumulator
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto cov_matrix_constrain_lkj(T&& x, size_t k, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [&lp, k](auto&& v) {
    return cov_matrix_constrain_lkj<Jacobian>(std::forward<decltype(v)>(v), k,
                                              lp);
  });
}

}  // namespace math
}  // namespace stan

#endif
