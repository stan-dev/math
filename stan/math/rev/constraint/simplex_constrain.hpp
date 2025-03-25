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
 * The transform is based on a centered stick-breaking process.
 *
 * @tparam T Type of vector to constrain
 * @param y Free vector input of dimensionality K - 1
 * @return Simplex of dimensionality K
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto simplex_constrain(const T& y) {
  using ret_type = plain_type_t<T>;

  size_t N = y.size();
  arena_t<T> arena_y = y;

  arena_t<ret_type> arena_x = simplex_constrain(arena_y.val());

  if (unlikely(N == 0)) {
    return ret_type(arena_x);
  }

  reverse_pass_callback([arena_y, arena_x]() mutable {
    int N = arena_y.size();
    double stick_len_val = arena_x.val().coeff(N);
    double stick_len_adj = arena_x.adj().coeff(N);
    for (Eigen::Index k = N; k-- > 0;) {
      double z_k = inv_logit(arena_y.val().coeff(k) - log(N - k));
      arena_x.adj().coeffRef(k) -= stick_len_adj;
      stick_len_val += arena_x.val().coeff(k);
      stick_len_adj += arena_x.adj().coeff(k) * z_k;
      double arena_z_adj = arena_x.adj().coeff(k) * stick_len_val;
      arena_y.adj().coeffRef(k) += arena_z_adj * z_k * (1.0 - z_k);
    }
  });

  return ret_type(arena_x);
}

/**
 * Return the simplex corresponding to the specified free vector
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam T type of the vector to constrain
 * @param y Free vector input of dimensionality N.
 * @param lp Log probability reference to increment.
 * @return Simplex of dimensionality N + 1.
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
auto simplex_constrain(const T& y, scalar_type_t<T>& lp) {
  using ret_type = plain_type_t<T>;

  size_t N = y.size();
  arena_t<T> arena_y = y;

  arena_t<ret_type> arena_x = simplex_constrain(arena_y.val(), lp);

  if (unlikely(N == 0)) {
    return ret_type(arena_x);
  }

  reverse_pass_callback([arena_y, arena_x, lp]() mutable {
    int N = arena_y.size();
    double stick_len_val = arena_x.val().coeff(N);
    double stick_len_adj = arena_x.adj().coeff(N);
    for (Eigen::Index k = N; k-- > 0;) {
      arena_x.adj().coeffRef(k) -= stick_len_adj;
      stick_len_val += arena_x.val().coeff(k);
      double log_N_minus_k = std::log(N - k);
      double adj_y_k = arena_y.val().coeff(k) - log_N_minus_k;
      double z_k = inv_logit(adj_y_k);

      arena_y.adj().coeffRef(k) -= lp.adj() * z_k;
      arena_y.adj().coeffRef(k) += lp.adj() * inv_logit(-adj_y_k);
      stick_len_adj += lp.adj() / stick_len_val;
      stick_len_adj += arena_x.adj().coeff(k) * z_k;
      double arena_z_adj = arena_x.adj().coeff(k) * stick_len_val;
      arena_y.adj().coeffRef(k) += arena_z_adj * z_k * (1.0 - z_k);
    }
  });

  return ret_type(arena_x);
}

}  // namespace math
}  // namespace stan
#endif
