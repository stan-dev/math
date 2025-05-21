#ifndef STAN_MATH_REV_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/constraint/sum_to_zero_constrain.hpp>
#include <stan/math/prim/constraint/simplex_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the simplex corresponding to the specified free vector.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A vector with (K-1) unconstrained values
 * will produce a simplex of size K.
 *
 * The simplex transform is defined using the inverse of the
 * isometric log ratio (ILR) transform. This is equivalent to
 * `softmax(sum_to_zero_constrain(y))`.
 *
 * @tparam T Type of vector to constrain
 * @param y Free vector input of dimensionality K - 1
 * @return Simplex of dimensionality K
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto simplex_constrain(const T& y) {
  using ret_type = plain_type_t<T>;

  const auto N = y.size();
  arena_t<T> arena_y = y;

  arena_t<ret_type> arena_x = simplex_constrain(arena_y.val());

  if (unlikely(N == 0)) {
    return ret_type(arena_x);
  }

  reverse_pass_callback([arena_y, arena_x]() mutable {
    auto&& res_val = arena_x.val();

    // backprop for softmax
    Eigen::VectorXd x_pre_softmax_adj = -res_val * arena_x.adj().dot(res_val)
                                        + res_val.cwiseProduct(arena_x.adj());

    // backprop for sum_to_zero_constrain
    internal::sum_to_zero_vector_backprop(arena_y.adj(), x_pre_softmax_adj);
  });

  return ret_type(arena_x);
}

/**
 * Return the simplex corresponding to the specified free vector
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined using the inverse of the
 * isometric log ratio (ILR) transform. This is equivalent to
 * `softmax(sum_to_zero_constrain(y))`.
 *
 * @tparam T type of the vector to constrain
 * @param y Free vector input of dimensionality N.
 * @param lp Log probability reference to increment.
 * @return Simplex of dimensionality N + 1.
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto simplex_constrain(const T& y, scalar_type_t<T>& lp) {
  using ret_type = plain_type_t<T>;

  const auto N = y.size();
  arena_t<T> arena_y = y;

  double lp_val = 0.0;
  arena_t<ret_type> arena_x = simplex_constrain(arena_y.val(), lp_val);
  lp += lp_val;

  if (unlikely(N == 0)) {
    return ret_type(arena_x);
  }

  reverse_pass_callback([arena_y, arena_x, lp]() mutable {
    auto&& res_val = arena_x.val();

    // backprop for log jacobian contribution to log density is equivalent to
    // arena_x.adj().array() += lp.adj() / res_val.array();
    // but is folded into the following to avoid needing to modify the adjoints
    // in-place

    // backprop for softmax
    Eigen::VectorXd x_pre_softmax_adj
        = -res_val * (arena_x.adj().dot(res_val) + res_val.size() * lp.adj())
          + (res_val.cwiseProduct(arena_x.adj()).array() + lp.adj()).matrix();

    // backprop for sum_to_zero_constrain
    internal::sum_to_zero_vector_backprop(arena_y.adj(), x_pre_softmax_adj);
  });

  return ret_type(arena_x);
}

}  // namespace math
}  // namespace stan
#endif
