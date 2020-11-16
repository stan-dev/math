#ifndef STAN_MATH_REV_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an increasing positive ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @param x Free vector of scalars
 * @return Positive, increasing ordered vector
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto positive_ordered_constrain(const T& x) {
  using ret_type = plain_type_t<T>;

  size_t N = x.size();

  if (unlikely(N == 0)) {
    return ret_type(x);
  }

  arena_t<T> arena_x = x;
  Eigen::VectorXd y_val(N);
  arena_t<Eigen::VectorXd> exp_x(N);

  exp_x(0) = std::exp(value_of(arena_x)(0));
  y_val(0) = exp_x.coeff(0);
  for (int n = 1; n < N; ++n) {
    exp_x(n) = std::exp(value_of(arena_x)(n));
    y_val(n) = y_val.coeff(n - 1) + exp_x.coeff(n);
  }

  arena_t<ret_type> y = y_val;

  reverse_pass_callback([arena_x, exp_x, y]() mutable {
    const size_t N = arena_x.size();
    double rolling_adjoint_sum = 0.0;

    for (int n = N; --n >= 0;) {
      rolling_adjoint_sum += y.adj().coeff(n);
      arena_x.adj().coeffRef(n) += exp_x.coeff(n) * rolling_adjoint_sum;
    }
  });

  return ret_type(y);
}

}  // namespace math
}  // namespace stan
#endif
