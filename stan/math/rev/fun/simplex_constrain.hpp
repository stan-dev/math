#ifndef STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
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
template<typename T,
	 require_eigen_col_vector_vt<is_var, T>* = nullptr>
inline auto simplex_constrain(const T& y) {
  using ret_type = plain_type_t<T>;
  
  size_t N = y.size();

  arena_t<T> arena_y = y;
  
  arena_t<Eigen::VectorXd> diag(N);
  arena_t<Eigen::VectorXd> z(N);

  Eigen::VectorXd x_val(N + 1);

  double stick_len(1.0);
  for (int k = 0; k < N; ++k) {
    double log_N_minus_k = std::log(N - k);
    z.coeffRef(k) = inv_logit(value_of(arena_y.coeff(k)) - log_N_minus_k);
    diag.coeffRef(k) = stick_len * z.coeff(k)
                       * inv_logit(log_N_minus_k - value_of(arena_y.coeff(k)));
    x_val.coeffRef(k) = stick_len * z.coeff(k);
    stick_len -= x_val(k);
  }
  x_val.coeffRef(N) = stick_len;

  arena_t<ret_type> x = x_val;

  if(N > 0) {
    reverse_pass_callback([arena_y, x, diag, z]() mutable {
      size_t N = arena_y.size();
      double acc = x.adj().coeff(N);

      arena_y.adj().coeffRef(N - 1)
	+= diag.coeff(N - 1) * (x.adj().coeff(N - 1) - acc);
      for (int n = N - 1; --n >= 0;) {
	acc = x.adj().coeff(n + 1) * z.coeff(n + 1)
	  + (1 - z.coeff(n + 1)) * acc;
	arena_y.adj().coeffRef(n) += diag.coeff(n) * (x.adj().coeff(n) - acc);
      }
    });
  }

  return ret_type(x);
}

}  // namespace math
}  // namespace stan
#endif
