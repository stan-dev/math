#ifndef STAN_MATH_REV_FUN_STOCHASTIC_ROW_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_STOCHASTIC_ROW_CONSTRAIN_HPP

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
 * Return a row stochastic matrix.
 *
 * @tparam T Type of matrix to constrain
 * @param y Free vector input of dimensionality (N, K - 1)
 * @return Matrix with Simplexes along the rows of dimensionality (N, K)
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline plain_type_t<T> stochastic_row_constrain(const T& y) {
  using ret_type = plain_type_t<T>;
  const Eigen::Index N = y.rows();
  const Eigen::Index M = y.cols();
  arena_t<Eigen::MatrixXd> x_val(N, M + 1);
  if (unlikely(N == 0 || M == 0)) {
    return ret_type(x_val);
  }
  arena_t<T> arena_y = y;
  arena_t<Eigen::MatrixXd> arena_z(N, M);
  Eigen::Array<double, -1, 1> stick_len = Eigen::Array<double, -1, 1>::Ones(N);
  for (Eigen::Index j = 0; j < M; ++j) {
    double log_N_minus_k = std::log(M - j);
    arena_z.col(j).array()
        = inv_logit((arena_y.col(j).val_op().array() - log_N_minus_k).matrix());
    x_val.col(j).array() = stick_len * arena_z.col(j).array();
    stick_len -= x_val.col(j).array();
  }
  x_val.col(M).array() = stick_len;
  arena_t<ret_type> arena_x = x_val;
  reverse_pass_callback([arena_y, arena_x, arena_z]() mutable {
    const Eigen::Index M = arena_y.cols();
    auto arena_y_arr = arena_y.array();
    auto arena_x_arr = arena_x.array();
    auto arena_z_arr = arena_z.array();
    auto stick_len_val_arr = arena_x_arr.col(M).val_op().eval();
    auto stick_len_adj_arr = arena_x_arr.col(M).adj_op().eval();
    for (Eigen::Index k = M; k-- > 0;) {
      arena_x_arr.col(k).adj() -= stick_len_adj_arr;
      stick_len_val_arr += arena_x_arr.col(k).val_op();
      stick_len_adj_arr += arena_x_arr.col(k).adj_op() * arena_z_arr.col(k);
      arena_y_arr.col(k).adj() += arena_x_arr.adj_op().col(k)
                                  * stick_len_val_arr * arena_z_arr.col(k)
                                  * (1.0 - arena_z_arr.col(k));
    }
  });
  return ret_type(arena_x);
}

/**
 * Return a row stochastic matrix
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam T type of the matrix to constrain
 * @param y Free matrix input of dimensionality (N, K).
 * @param lp Log probability reference to increment.
 * @return Matrix with simplexes along the rows of dimensionality (N, K + 1).
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline plain_type_t<T> stochastic_row_constrain(const T& y,
                                                scalar_type_t<T>& lp) {
  using ret_type = plain_type_t<T>;
  const Eigen::Index N = y.rows();
  const Eigen::Index M = y.cols();
  arena_t<Eigen::MatrixXd> x_val(N, M + 1);
  if (unlikely(N == 0 || M == 0)) {
    return ret_type(x_val);
  }
  arena_t<T> arena_y = y;
  arena_t<Eigen::MatrixXd> arena_z(N, M);
  Eigen::Array<double, -1, 1> stick_len = Eigen::Array<double, -1, 1>::Ones(N);
  for (Eigen::Index j = 0; j < M; ++j) {
    double log_N_minus_k = std::log(M - j);
    auto adj_y_k = arena_y.col(j).val_op().array() - log_N_minus_k;
    arena_z.col(j).array() = inv_logit(adj_y_k);
    x_val.col(j).array() = stick_len * arena_z.col(j).array();
    lp += sum(log(stick_len)) - sum(log1p_exp(-adj_y_k))
          - sum(log1p_exp(adj_y_k));
    stick_len -= x_val.col(j).array();
  }
  x_val.col(M).array() = stick_len;
  arena_t<ret_type> arena_x = x_val;
  reverse_pass_callback([arena_y, arena_x, arena_z, lp]() mutable {
    const Eigen::Index M = arena_y.cols();
    auto arena_y_arr = arena_y.array();
    auto arena_x_arr = arena_x.array();
    auto arena_z_arr = arena_z.array();
    auto stick_len_val = arena_x_arr.col(M).val_op().eval();
    auto stick_len_adj = arena_x_arr.col(M).adj_op().eval();
    for (Eigen::Index k = M; k-- > 0;) {
      const double log_N_minus_k = std::log(M - k);
      arena_x_arr.col(k).adj() -= stick_len_adj;
      stick_len_val += arena_x_arr.col(k).val_op();
      stick_len_adj += lp.adj() / stick_len_val
                       + arena_x_arr.adj_op().col(k) * arena_z_arr.col(k);
      auto adj_y_k = arena_y_arr.col(k).val_op() - log_N_minus_k;
      arena_y_arr.col(k).adj()
          += -(lp.adj() * inv_logit(adj_y_k)) + lp.adj() * inv_logit(-adj_y_k)
             + arena_x_arr.col(k).adj_op() * stick_len_val * arena_z_arr.col(k)
                   * (1.0 - arena_z_arr.col(k));
    }
  });
  return ret_type(arena_x);
}

}  // namespace math
}  // namespace stan
#endif
