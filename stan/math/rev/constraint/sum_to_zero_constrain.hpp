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

    double sum_u_adj = 0;
    for (int i = 0; i < N; ++i) {
      double n = i + 1;

      // adjoint of the reverse cumulative sum computed in the forward mode
      sum_u_adj += arena_z.adj()(i);

      // adjoint of the offset subtraction
      double v_adj = -arena_z.adj()(i + 1) * n;

      double w_adj = v_adj + sum_u_adj;

      arena_y.adj()(i) += w_adj / sqrt(n * (n + 1));
    }
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
