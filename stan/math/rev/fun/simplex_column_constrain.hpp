#ifndef STAN_MATH_REV_FUN_SIMPLEX_COLUMN_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_SIMPLEX_COLUMN_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the simplex corresponding to the specified free vector.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A vector with (K-1) unconstrained values
 * will produce a simplex of size K.
 *
 * The transform is based on a centered stick-breaking process.
 *
 * @tparam T Type of vector to constrain
 * @param y Free vector input of dimensionality K - 1
 * @return Simplex of dimensionality K
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto simplex_column_constrain(const T& y) {
  using ret_type = plain_type_t<T>;

  const Eigen::Index N = y.rows();
  const Eigen::Index M = y.cols();
  arena_t<T> arena_y = y;
  arena_t<Eigen::MatrixXd> arena_z(N, M);
  arena_t<Eigen::MatrixXd> x_val(N + 1, M);
  for (Eigen::Index j = 0; j < x_val.cols(); ++j) {
    double stick_len(1.0);
    for (Eigen::Index k = 0; k < N; ++k) {
      double log_N_minus_k = std::log(N - k);
      arena_z.coeffRef(k, j) = inv_logit(arena_y.val().coeff(k, j) - log_N_minus_k);
      x_val.coeffRef(k, j) = stick_len * arena_z.coeff(k, j);
      stick_len -= x_val.coeff(k, j);
    }
    x_val.coeffRef(N, j) = stick_len;
  }

  arena_t<ret_type> arena_x = x_val;

  if (unlikely(N == 0)) {
    return ret_type(arena_x);
  }

  reverse_pass_callback([arena_y, arena_x, arena_z]() mutable {
    const Eigen::Index N = arena_y.rows();
    const Eigen::Index M = arena_y.cols();
    for (Eigen::Index j = 0; j < M; ++j) {
      double stick_len_val = arena_x.val().coeff(N, j);
      double stick_len_adj = arena_x.adj().coeff(N, j);
      for (Eigen::Index k = N; k-- > 0;) {
        arena_x.adj().coeffRef(k, j) -= stick_len_adj;
        stick_len_val += arena_x.val().coeff(k, j);
        stick_len_adj += arena_x.adj().coeff(k, j) * arena_z.coeff(k, j);
        double arena_z_adj = arena_x.adj().coeff(k, j) * stick_len_val;
        arena_y.adj().coeffRef(k, j)
            += arena_z_adj * arena_z.coeff(k, j) * (1.0 - arena_z.coeff(k, j));
      }
    }
  });

  return ret_type(arena_x);
}

/**
 * Return the simplex corresponding to the specified free vector
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam T type of the vector to constrain
 * @param y Free vector input of dimensionality N.
 * @param lp Log probability reference to increment.
 * @return Simplex of dimensionality N + 1.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
auto simplex_column_constrain(const T& y, scalar_type_t<T>& lp) {
  using ret_type = plain_type_t<T>;

  const Eigen::Index N = y.rows();
  const Eigen::Index M = y.cols();
  arena_t<T> arena_y = y;
  arena_t<Eigen::MatrixXd> arena_z(N, M);
  arena_t<Eigen::MatrixXd> x_val(N + 1, M);
  for (Eigen::Index j = 0; j < M; ++j) {
    double stick_len(1.0);
    for (Eigen::Index k = 0; k < N; ++k) {
      double log_N_minus_k = std::log(N - k);
      double adj_y_k = arena_y.val().coeff(k, j) - log_N_minus_k;
      arena_z.coeffRef(k, j) = inv_logit(adj_y_k);
      x_val.coeffRef(k, j) = stick_len * arena_z.coeff(k, j);
      lp += log(stick_len);
      lp -= log1p_exp(-adj_y_k);
      lp -= log1p_exp(adj_y_k);
      stick_len -= x_val(k, j);
    }
    x_val.coeffRef(N, j) = stick_len;
  }

  arena_t<ret_type> arena_x = x_val;

  if (unlikely(N == 0)) {
    return ret_type(arena_x);
  }

  reverse_pass_callback([arena_y, arena_x, arena_z, lp]() mutable {
    const Eigen::Index N = arena_y.rows();
    const Eigen::Index M = arena_y.cols();
    for (Eigen::Index j = 0; j < M; ++j) {
      double stick_len_val = arena_x.val().coeff(N, j);
      double stick_len_adj = arena_x.adj().coeff(N, j);
      for (Eigen::Index k = N; k-- > 0;) {
        arena_x.adj().coeffRef(k, j) -= stick_len_adj;
        stick_len_val += arena_x.val().coeff(k, j);
        double log_N_minus_k = std::log(N - k);
        double adj_y_k = arena_y.val().coeff(k, j) - log_N_minus_k;
        arena_y.adj().coeffRef(k, j) -= lp.adj() * inv_logit(adj_y_k);
        arena_y.adj().coeffRef(k, j) += lp.adj() * inv_logit(-adj_y_k);
        stick_len_adj += lp.adj() / stick_len_val;
        stick_len_adj += arena_x.adj().coeff(k, j) * arena_z.coeff(k, j);
        double arena_z_adj = arena_x.adj().coeff(k, j) * stick_len_val;
        arena_y.adj().coeffRef(k, j)
            += arena_z_adj * arena_z.coeff(k, j) * (1.0 - arena_z.coeff(k, j));
      }

    }
  });

  return ret_type(arena_x);
}

}  // namespace math
}  // namespace stan
#endif
