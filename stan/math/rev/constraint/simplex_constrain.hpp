#ifndef STAN_MATH_REV_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_SIMPLEX_CONSTRAIN_HPP

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
 * The simplex transform is defined using the Isometric Log Ratio
 * transform (ILR). This is equilvalent to
 * `softmax(sum_to_zero_constrain(y))`.
 *
 * @tparam T Type of vector to constrain
 * @param y Free vector input of dimensionality K - 1
 * @return Simplex of dimensionality K
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto simplex_constrain(const T& y) {
  return softmax(sum_to_zero_constrain(y));
}

/**
 * Return the simplex corresponding to the specified free vector
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined using the Isometric Log Ratio
 * transform (ILR). This is equilvalent to
 * `softmax(sum_to_zero_constrain(y))`.
 *
 * @tparam T type of the vector to constrain
 * @param y Free vector input of dimensionality N.
 * @param lp Log probability reference to increment.
 * @return Simplex of dimensionality N + 1.
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
auto simplex_constrain(const T& y, scalar_type_t<T>& lp) {
  using ret_type = plain_type_t<T>;

  arena_t<ret_type> log_x = log_softmax(sum_to_zero_constrain(y));

  const auto N = y.size() + 1;

  lp += sum(log_x) + 0.5 * log(N);

  return ret_type(exp(log_x));
}

}  // namespace math
}  // namespace stan
#endif
