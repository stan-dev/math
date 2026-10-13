#ifndef STAN_MATH_REV_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_CHOLESKY_CORR_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/log_cosh.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <stan/math/prim/constraint/cholesky_corr_constrain.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

namespace internal {
/**
 * Return the Cholesky factor of a correlation matrix from unconstrained
 * values `y`, whose partial correlations are `tanh(y)`. See the `prim`
 * version for the formulas.
 *
 * Reverse pass: in row `r`, element `c` depends on its own `y` through
 * `tanh(y) exp(-a)`, with derivative `exp(-2 log(cosh(y)) - a)`, and every
 * later element of the row, the diagonal included, depends on it through
 * `exp(-a)`, with derivative `-tanh(y)` times the element.
 *
 * @tparam ColumnOrder order of the partial correlations in `y`, see
 * `cholesky_corr_position`
 * @tparam T type of the `var_value` vector
 * @param y unconstrained values, `K (K - 1) / 2` of them
 * @param K number of rows
 * @return Cholesky factor of a correlation matrix
 */
template <bool ColumnOrder, typename T, require_var_vector_t<T>* = nullptr>
inline var_value<Eigen::MatrixXd> cholesky_corr_factor(const T& y,
                                                       Eigen::Index K) {
  using Eigen::Index;
  if (K == 0) {
    return Eigen::MatrixXd();
  }
  const Index n = y.size();
  arena_t<Eigen::VectorXd> tanh_y(n);
  // derivative of each element with respect to its own unconstrained value
  arena_t<Eigen::VectorXd> d_own(n);
  arena_t<Eigen::MatrixXd> L_val = Eigen::MatrixXd::Zero(K, K);
  L_val.coeffRef(0, 0) = 1;
  for (Index r = 1; r < K; ++r) {
    // exp(-a), as the product of sech(y) over the earlier elements of the row
    double exp_minus_a = 1.0;
    for (Index c = 0; c < r; ++c) {
      const Index k = cholesky_corr_position<ColumnOrder>(r, c, K);
      const std::pair<double, double> t_s = tanh_and_sech(y.val().coeff(k));
      tanh_y.coeffRef(k) = t_s.first;
      d_own.coeffRef(k) = t_s.second * t_s.second * exp_minus_a;
      L_val.coeffRef(r, c) = t_s.first * exp_minus_a;
      exp_minus_a *= t_s.second;
    }
    L_val.coeffRef(r, r) = exp_minus_a;
  }
  var_value<Eigen::MatrixXd> L = L_val;
  reverse_pass_callback([y, L, tanh_y, d_own, K]() mutable {
    for (Index r = 1; r < K; ++r) {
      double later = L.adj().coeff(r, r) * L.val().coeff(r, r);
      for (Index c = r - 1; c >= 0; --c) {
        const Index k = cholesky_corr_position<ColumnOrder>(r, c, K);
        const double adj_rc = L.adj().coeff(r, c);
        y.adj().coeffRef(k)
            += adj_rc * d_own.coeff(k) - tanh_y.coeff(k) * later;
        later += adj_rc * L.val().coeff(r, c);
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
inline var_value<Eigen::MatrixXd> cholesky_corr_constrain(const T& y, int K) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  return internal::cholesky_corr_factor<false>(y, K);
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
    const T& y, int K, scalar_type_t<T>& lp) {
  int k_choose_2 = (K * (K - 1)) / 2;
  check_size_match("cholesky_corr_constrain", "y.size()", y.size(),
                   "k_choose_2", k_choose_2);
  if (K == 0) {
    return Eigen::MatrixXd();
  }
  lp -= dot_product(internal::cholesky_corr_log_jacobian_coefficients(K),
                    log_cosh(y));
  return internal::cholesky_corr_factor<false>(y, K);
}

}  // namespace math
}  // namespace stan
#endif
