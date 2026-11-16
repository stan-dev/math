#ifndef STAN_MATH_REV_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/log_cosh.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <stan/math/prim/fun/sech.hpp>
#include <stan/math/prim/fun/tanh.hpp>
#include <stan/math/prim/constraint/cholesky_corr_constrain.hpp>
#include <utility>

namespace stan {
namespace math {

namespace internal {
/**
 * Return the Cholesky factor of a correlation matrix from unconstrained
 * values `y`, whose partial correlations are `tanh(y)`, as the `prim`
 * version does. The reverse pass propagates each adjoint to its own `y` and,
 * through the shared `exp(-a)`, to every earlier `y` of the row.
 *
 * @tparam Order order of the partial correlations in `y`
 * @tparam T type of the `var_value` vector
 * @param y unconstrained values, `K (K - 1) / 2` of them
 * @param K number of rows
 * @return Cholesky factor of a correlation matrix
 */
template <CorrOrder Order, typename T, require_var_vector_t<T>* = nullptr>
inline var_value<Eigen::MatrixXd> cholesky_corr_factor(T&& y, Eigen::Index K) {
  using Eigen::Index;
  if (K == 0) {
    return Eigen::MatrixXd();
  }
  auto arena_y = to_arena(std::forward<T>(y));
  arena_t<Eigen::VectorXd> tanh_y(arena_y.size());
  // derivative of each element with respect to its own unconstrained value
  arena_t<Eigen::VectorXd> d_own(arena_y.size());
  arena_t<Eigen::MatrixXd> L_val = Eigen::MatrixXd::Zero(K, K);
  L_val.coeffRef(0, 0) = 1;
  for (Index row = 1; row < K; ++row) {
    // exp(-acc), as the product of sech(y) over the earlier elements of the
    // row
    double exp_neg_acc = 1.0;
    for (Index col = 0; col < row; ++col) {
      const Index pos = cholesky_corr_position<Order>(row, col, K);
      const double y_val = arena_y.val().coeff(pos);
      const double sech_val = sech(y_val);
      tanh_y.coeffRef(pos) = tanh(y_val);
      d_own.coeffRef(pos) = sech_val * sech_val * exp_neg_acc;
      L_val.coeffRef(row, col) = tanh_y.coeff(pos) * exp_neg_acc;
      exp_neg_acc *= sech_val;
    }
    L_val.coeffRef(row, row) = exp_neg_acc;
  }
  var_value<Eigen::MatrixXd> L = L_val;
  reverse_pass_callback([arena_y, L, tanh_y, d_own, K]() mutable {
    for (Index row = 1; row < K; ++row) {
      // sum of adj * val over the elements of the row after column col, which
      // all depend on y_col through the shared exp(-acc)
      double running_sum = L.adj().coeff(row, row) * L.val().coeff(row, row);
      for (Index col = row - 1; col >= 0; --col) {
        const Index pos = cholesky_corr_position<Order>(row, col, K);
        const double adj_rc = L.adj().coeff(row, col);
        arena_y.adj().coeffRef(pos)
            += adj_rc * d_own.coeff(pos) - tanh_y.coeff(pos) * running_sum;
        running_sum += adj_rc * L.val().coeff(row, col);
      }
    }
  });
  return L;
}
}  // namespace internal

/**
 * Return the Cholesky factor of the correlation matrix of the specified
 * size read from the unconstrained vector `y`. A total of K choose 2
 * elements are required to build a K by K Cholesky factor.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param y Vector of unconstrained values
 * @param K number of rows
 * @return Cholesky factor of correlation matrix
 */
template <typename T, require_var_vector_t<T>* = nullptr>
inline var_value<Eigen::MatrixXd> cholesky_corr_constrain(T&& y, int K) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  return internal::cholesky_corr_factor<internal::CorrOrder::row>(
      std::forward<T>(y), K);
}

/**
 * Return the Cholesky factor of the correlation matrix of the specified
 * size read from the unconstrained vector `y`. A total of K choose 2
 * elements are required to build a K by K Cholesky factor.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param y Vector of unconstrained values
 * @param K number of rows
 * @param[out] lp Log density that is incremented with the log Jacobian
 * @return Cholesky factor of correlation matrix
 */
template <typename T, require_var_vector_t<T>* = nullptr>
inline var_value<Eigen::MatrixXd> cholesky_corr_constrain(
    T&& y, int K, scalar_type_t<T>& lp) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  lp -= dot_product(internal::cholesky_corr_log_jacobian_coefficients(K),
                    log_cosh(y));
  return internal::cholesky_corr_factor<internal::CorrOrder::row>(
      std::forward<T>(y), K);
}

}  // namespace math
}  // namespace stan
#endif
