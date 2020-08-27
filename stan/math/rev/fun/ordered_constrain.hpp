#ifndef STAN_MATH_REV_FUN_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_ORDERED_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an increasing ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @param x Free vector of scalars
 * @return Increasing ordered vector
 */
inline Eigen::Matrix<var, Eigen::Dynamic, 1> ordered_constrain(
    const Eigen::Matrix<var, Eigen::Dynamic, 1>& x) {
  using std::exp;

  size_t N = x.size();
  if (N == 0) {
    return x;
  }

  Eigen::VectorXd y_val(N);
  arena_matrix<Eigen::VectorXd> exp_x(N - 1);

  y_val.coeffRef(0) = value_of(x.coeff(0));
  for (int n = 1; n < N; ++n) {
    exp_x.coeffRef(n - 1) = exp(value_of(x.coeff(n)));
    y_val.coeffRef(n) = y_val.coeff(n - 1) + exp_x.coeff(n - 1);
  }

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> y = y_val;
  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> arena_x = x;

  reverse_pass_callback([=]() mutable {
    double rolling_adjoint_sum = 0.0;

    if (N > 0) {
      for (int n = N - 1; n > 0; --n) {
        rolling_adjoint_sum += y.adj().coeff(n);
        arena_x.adj().coeffRef(n) += exp_x.coeff(n - 1) * rolling_adjoint_sum;
      }
      arena_x.adj().coeffRef(0) += rolling_adjoint_sum + y.adj().coeff(0);
    }
  });

  return y;
}

}  // namespace math
}  // namespace stan
#endif
