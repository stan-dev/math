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

  size_t N = y.size();
  Eigen::VectorXd x_val(N + 1);

  arena_t<T> arena_y = y;

  if (unlikely(N == 0)) {
    x_val << 0;
    return ret_type(x_val);
  }

  x_val.head(N) = y.val();
  arena_t<ret_type> arena_x = x_val;

  var x_N = -sum(y);

  arena_x.coeffRef(N) = x_N;

  reverse_pass_callback([arena_y, arena_x, x_N, N]() mutable {
    arena_y.adj() += arena_x.adj().head(N);
    x_N.adj() += arena_x.adj().coeff(N);
  });

  return ret_type(arena_x);
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
auto sum_to_zero_constrain(const T& y, scalar_type_t<T>& lp) {
  return sum_to_zero_constrain(y);
}

}  // namespace math
}  // namespace stan
#endif
