#ifndef STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/constraint/sum_to_zero_constrain.hpp>
#include <stan/math/prim/fun/linspaced_vector.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined such that the first K-1
 * elements are unconstrained and the last element is the negative
 * sum of those elements.
 *
 * @tparam T type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Zero-sum vector of dimensionality K.
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto sum_to_zero_constrain(const T& y) {
  using ret_type = plain_type_t<T>;
  if (unlikely(y.size() == 0)) {
    return arena_t<ret_type>(Eigen::VectorXd{{0}});
  }
  auto arena_y = to_arena(y);
  arena_t<ret_type> arena_z = sum_to_zero_constrain(arena_y.val());

  reverse_pass_callback([arena_y, arena_z]() mutable {
    const auto N = arena_y.size();

    Eigen::VectorXd ns = linspaced_vector(N, 1, N );

    /*
    u3.adj += z.adj[1:N-1]
    v.adj -= z.adj[2:N]
    w.adj += v.adj .* ns
    u2.adj += reverse(u3.adj)
    u1.adj += cumulative_sum(u2.adj)
    w.adj += reverse(u1.adj)
    y.adj += w.adj ./ sqrt(ns .* (ns + 1))
    */

    Eigen::VectorXd u_adj = arena_z.adj_op().head(N).eval();
    Eigen::VectorXd v_adj = -arena_z.adj_op().tail(N).eval();

    Eigen::VectorXd w_from_v = v_adj.array() * ns.array();

    Eigen::VectorXd w_from_u = cumulative_sum(u_adj);

    Eigen::VectorXd w_adj = (w_from_v.array() + w_from_u.array());

    arena_y.adj() += (w_adj.array() / sqrt(ns.array() * (ns.array() + 1))).matrix();
  });

  return arena_z;
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined such that the first K-1
 * elements are unconstrained and the last element is the negative
 * sum of those elements. This is a linear transform, with no
 * Jacobian.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @param lp unused
 * @return Zero-sum vector of dimensionality K.
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto sum_to_zero_constrain(const T& y, scalar_type_t<T>& lp) {
  return sum_to_zero_constrain(y);
}

}  // namespace math
}  // namespace stan
#endif
