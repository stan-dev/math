#ifndef STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
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
 * The transform is based on a centered stick-breaking process.
 *
 * @param y Free vector input of dimensionality K - 1
 * @return Simplex of dimensionality K
 */
inline Eigen::Matrix<var, Eigen::Dynamic, 1> simplex_constrain(
    const Eigen::Matrix<var, Eigen::Dynamic, 1>& y) {
  size_t N = y.size();

  arena_matrix<Eigen::VectorXd> diag(N);
  arena_matrix<Eigen::VectorXd> z(N);

  Eigen::VectorXd x_val(N + 1);

  double stick_len(1.0);
  for (int k = 0; k < N; ++k) {
    double log_N_minus_k = std::log(N - k);
    z.coeffRef(k) = inv_logit(value_of(y.coeffRef(k)) - log_N_minus_k);
    diag.coeffRef(k) = stick_len * z.coeffRef(k) * inv_logit(log_N_minus_k - value_of(y.coeffRef(k)));
    x_val.coeffRef(k) = stick_len * z.coeffRef(k);
    stick_len -= x_val(k);
  }
  x_val.coeffRef(N) = stick_len;

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> x = x_val;
  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> arena_y = y;

  reverse_pass_callback([=]() mutable {
    const auto& adj = x.adj();
    double acc = adj(N);

    if (N > 0) {
      arena_y.adj().coeffRef(N - 1) += diag.coeffRef(N - 1) * (adj(N - 1) - acc);
      for (int n = N - 1; --n >= 0;) {
	acc = adj(n + 1) * z.coeffRef(n + 1) + (1 - z.coeffRef(n + 1)) * acc;
	arena_y.adj().coeffRef(n) += diag.coeffRef(n) * (adj(n) - acc);
      }
    }
  });

  return x;
}

}  // namespace math
}  // namespace stan
#endif
