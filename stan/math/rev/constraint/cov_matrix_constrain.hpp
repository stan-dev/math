#ifndef STAN_MATH_REV_CONSTRAINT_COV_MATRIX_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_COV_MATRIX_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/rev/fun/multiply_lower_tri_self_transpose.hpp>
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
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x The vector to convert to a covariance matrix.
 * @param K The number of rows and columns of the resulting
 * covariance matrix.
 * @throws std::invalid_argument if (x.size() != K + (K choose 2)).
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> cov_matrix_constrain(const T& x, Eigen::Index K) {
  using std::exp;

  check_size_match("cov_matrix_constrain", "x.size()", x.size(),
                   "K + (K choose 2)", (K * (K + 1)) / 2);
  arena_t<Eigen::MatrixXd> L_val = Eigen::MatrixXd::Zero(K, K);
  int i = 0;
  for (Eigen::Index m = 0; m < K; ++m) {
    L_val.row(m).head(m) = x.val().segment(i, m);
    i += m;
    L_val.coeffRef(m, m) = exp(x.val().coeff(i));
    i++;
  }

  var_value<Eigen::MatrixXd> L = L_val;

  reverse_pass_callback([x, L]() mutable {
    Eigen::Index K = L.rows();
    int i = x.size();
    for (int m = K - 1; m >= 0; --m) {
      i--;
      x.adj().coeffRef(i) += L.adj().coeff(m, m) * L.val().coeff(m, m);
      i -= m;
      x.adj().segment(i, m) += L.adj().row(m).head(m);
    }
  });

  return multiply_lower_tri_self_transpose(L);
}

/**
 * Return the symmetric, positive-definite matrix of dimensions K
 * by K resulting from transforming the specified finite vector of
 * size K plus (K choose 2).
 *
 * <p>See <code>cov_matrix_free()</code> for the inverse transform.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x The vector to convert to a covariance matrix.
 * @param K The dimensions of the resulting covariance matrix.
 * @param lp Reference
 * @throws std::domain_error if (x.size() != K + (K choose 2)).
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> cov_matrix_constrain(const T& x, Eigen::Index K,
                                                scalar_type_t<T>& lp) {
  using std::exp;
  using std::log;

  check_size_match("cov_matrix_constrain", "x.size()", x.size(),
                   "K + (K choose 2)", (K * (K + 1)) / 2);
  arena_t<Eigen::MatrixXd> L_val = Eigen::MatrixXd::Zero(K, K);
  int pos = 0;
  for (Eigen::Index m = 0; m < K; ++m) {
    L_val.row(m).head(m) = x.val().segment(pos, m);
    pos += m;
    L_val.coeffRef(m, m) = exp(x.val().coeff(pos));
    pos++;
  }

  // Jacobian for complete transform, including exp() above
  double lp_val = (K * LOG_TWO);
  for (Eigen::Index k = 0; k < K; ++k) {
    // only +1 because index from 0
    lp_val += (K - k + 1) * log(L_val.coeff(k, k));
  }

  lp += lp_val;

  var_value<Eigen::MatrixXd> L = L_val;

  reverse_pass_callback([x, L, lp]() mutable {
    Eigen::Index K = L.rows();
    for (Eigen::Index k = 0; k < K; ++k) {
      L.adj().coeffRef(k, k) += (K - k + 1) * lp.adj() / L.val().coeff(k, k);
    }
    int pos = x.size();
    for (int m = K - 1; m >= 0; --m) {
      pos--;
      x.adj().coeffRef(pos) += L.adj().coeff(m, m) * L.val().coeff(m, m);
      pos -= m;
      x.adj().segment(pos, m) += L.adj().row(m).head(m);
    }
  });

  return multiply_lower_tri_self_transpose(L);
}

}  // namespace math
}  // namespace stan

#endif
