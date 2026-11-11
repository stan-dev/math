#ifndef STAN_MATH_REV_CONSTRAINT_COV_MATRIX_CONSTRAIN_LKJ_HPP
#define STAN_MATH_REV_CONSTRAINT_COV_MATRIX_CONSTRAIN_LKJ_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/constraint/cholesky_corr_constrain.hpp>
#include <stan/math/rev/fun/diag_pre_multiply.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/log_cosh.hpp>
#include <stan/math/rev/fun/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/positive_constrain.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the `k x k` covariance matrix `(diag(sds) L) (diag(sds) L)^T`, with
 * `L` the Cholesky factor of the correlation matrix of the first
 * `k choose 2` unconstrained values and `sds = exp(x.tail(k))`. See
 * `corr_matrix_constrain` for the correlation transform.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations.
 * @param k Dimensionality of returned covariance matrix.
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <typename T, require_var_vector_t<T>* = nullptr>
inline var_value<Eigen::MatrixXd> cov_matrix_constrain_lkj(T&& x, size_t k) {
  size_t k_choose_2 = (k * (k - 1)) / 2;
  check_size_match("cov_matrix_constrain_lkj", "x.size()", x.size(),
                   "k_choose_2 + k", k_choose_2 + k);
  auto arena_x = to_arena(std::forward<T>(x));
  return multiply_lower_tri_self_transpose(diag_pre_multiply(
      positive_constrain(arena_x.tail(k)),
      internal::cholesky_corr_factor<internal::CorrOrder::column>(
          arena_x.head(k_choose_2), k)));
}

/**
 * Return the covariance matrix as `cov_matrix_constrain_lkj(x, k)` and add
 * the log absolute Jacobian determinant to `lp`. The Jacobian is the one of
 * `corr_matrix_constrain` plus `k log 2 + (k + 1) sum(x.tail(k))`.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x Input vector of unconstrained partial correlations and
 * standard deviations.
 * @param k Dimensionality of returned covariance matrix.
 * @param lp Log probability reference to increment.
 * @return Covariance matrix derived from the unconstrained partial
 * correlations and deviations.
 */
template <typename T, require_var_vector_t<T>* = nullptr>
inline var_value<Eigen::MatrixXd> cov_matrix_constrain_lkj(
    T&& x, size_t k, scalar_type_t<T>& lp) {
  size_t k_choose_2 = (k * (k - 1)) / 2;
  check_size_match("cov_matrix_constrain_lkj", "x.size()", x.size(),
                   "k_choose_2 + k", k_choose_2 + k);
  auto arena_x = to_arena(std::forward<T>(x));
  auto&& cpcs = arena_x.head(k_choose_2);
  auto&& log_sds = arena_x.tail(k);
  lp -= dot_product(internal::corr_matrix_log_jacobian_coefficients(k),
                    log_cosh(cpcs));
  lp += (k + 1.0) * sum(log_sds) + k * LOG_TWO;
  return multiply_lower_tri_self_transpose(diag_pre_multiply(
      positive_constrain(log_sds),
      internal::cholesky_corr_factor<internal::CorrOrder::column>(cpcs, k)));
}

}  // namespace math
}  // namespace stan

#endif
