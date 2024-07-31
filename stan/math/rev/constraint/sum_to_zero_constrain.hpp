#ifndef STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
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

  const auto N = y.size();
  if (unlikely(N == 0)) {
    return arena_t<ret_type>(Eigen::VectorXd{{0}});
  }
  Eigen::VectorXd x_val = Eigen::VectorXd::Zero(N + 1);
  auto arena_y = to_arena(y);
  double x_sum = -sum(arena_y.val());
  x_val.head(N) = arena_y.val();
  x_val(N) = x_sum;
  arena_t<ret_type> arena_x = x_val;
  reverse_pass_callback([arena_y, arena_x, x_sum, N]() mutable {
    arena_y.adj().array() -= arena_x.adj_op()(N);
    arena_y.adj() += arena_x.adj_op().head(N);
  });
  return arena_x;
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
