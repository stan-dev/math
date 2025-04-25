#ifndef STAN_MATH_REV_CONSTRAINT_STOCHASTIC_ROW_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_STOCHASTIC_ROW_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/constraint/stochastic_row_constrain.hpp>
#include <stan/math/rev/constraint/sum_to_zero_constrain.hpp>
#include <cmath>

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
inline auto stochastic_row_constrain(const T& y) {
  using ret_type = plain_type_t<T>;

  const auto N = y.rows();
  const auto M = y.cols();
  arena_t<T> arena_y = y;

  arena_t<ret_type> arena_x = stochastic_row_constrain(arena_y.val_op());

  if (unlikely(N == 0 || M == 0)) {
    return arena_x;
  }

  reverse_pass_callback([arena_y, arena_x]() mutable {
    const auto N = arena_y.rows();

    auto&& x_val = arena_x.val_op();
    auto&& x_adj = arena_x.adj_op();

    Eigen::VectorXd x_pre_softmax_adj(x_val.cols());
    for (Eigen::Index i = 0; i < N; ++i) {
      // backprop for softmax
      x_pre_softmax_adj.noalias()
          = -x_val.row(i) * x_adj.row(i).dot(x_val.row(i))
            + x_val.row(i).cwiseProduct(x_adj.row(i));

      // backprop for sum_to_zero_constrain
      internal::sum_to_zero_vector_backprop(arena_y.row(i).adj(),
                                            x_pre_softmax_adj);
    }
  });

  return arena_x;
}

/**
 * Return a row stochastic matrix
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined using the inverse of the
 * isometric log ratio (ILR) transform
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

  const auto N = y.rows();
  const auto M = y.cols();
  arena_t<T> arena_y = y;

  double lp_val = 0;
  arena_t<ret_type> arena_x
      = stochastic_row_constrain(arena_y.val_op(), lp_val);
  lp += lp_val;

  if (unlikely(N == 0 || M == 0)) {
    return arena_x;
  }

  reverse_pass_callback([arena_y, arena_x, lp]() mutable {
    const auto N = arena_y.rows();

    auto&& x_val = arena_x.val_op();
    auto&& x_adj = arena_x.adj_op();

    const auto x_val_cols = x_val.cols();

    Eigen::VectorXd x_pre_softmax_adj(x_val.cols());
    for (Eigen::Index i = 0; i < N; ++i) {
      // backprop for softmax
      x_pre_softmax_adj.noalias()
          = -x_val.row(i)
                * (x_adj.row(i).dot(x_val.row(i)) + lp.adj() * x_val_cols)
            + (x_val.row(i).cwiseProduct(x_adj.row(i)).array() + lp.adj())
                  .matrix();

      // backprop for sum_to_zero_constrain
      internal::sum_to_zero_vector_backprop(arena_y.row(i).adj(),
                                            x_pre_softmax_adj);
    }
  });

  return arena_x;
}

}  // namespace math
}  // namespace stan
#endif
