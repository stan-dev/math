#ifndef STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
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
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto simplex_constrain(const T& y) {
  using ret_type = plain_type_t<T>;

  size_t N = y.size();

  arena_t<T> arena_y = y;

  arena_t<Eigen::VectorXd> arena_diag(N);
  arena_t<Eigen::VectorXd> arena_z(N);

  Eigen::VectorXd x_val(N + 1);

  double stick_len(1.0);
  for (Eigen::Index k = 0; k < N; ++k) {
    double log_N_minus_k = std::log(N - k);
    arena_z.coeffRef(k) = inv_logit(arena_y.val().coeff(k) - log_N_minus_k);
    arena_diag.coeffRef(k)
        = stick_len * arena_z.coeff(k)
          * inv_logit(log_N_minus_k - arena_y.val().coeff(k));
    x_val.coeffRef(k) = stick_len * arena_z.coeff(k);
    stick_len -= x_val(k);
  }
  x_val.coeffRef(N) = stick_len;

  arena_t<ret_type> arena_x = x_val;
  if (unlikely(N == 0)) {
    return ret_type(arena_x);
  }

  reverse_pass_callback([arena_y, arena_x, arena_diag, arena_z]() mutable {
    auto N = arena_y.size();
    double acc = arena_x.adj().coeff(N);

    arena_y.adj().coeffRef(N - 1)
        += arena_diag.coeff(N - 1) * (arena_x.adj().coeff(N - 1) - acc);
    for (Eigen::Index n = N - 1; --n >= 0;) {
      acc = arena_x.adj().coeff(n + 1) * arena_z.coeff(n + 1)
            + (1 - arena_z.coeff(n + 1)) * acc;
      arena_y.adj().coeffRef(n)
          += arena_diag.coeff(n) * (arena_x.adj().coeff(n) - acc);
    }
  });

  return ret_type(arena_x);
}

}  // namespace math
}  // namespace stan
#endif
