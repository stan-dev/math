#ifndef STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 **/
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
plain_type_t<T> unit_vector_constrain(const T& y) {
  check_vector("unit_vector", "y", y);
  check_nonzero_size("unit_vector", "y", y);

  const auto& y_ref = to_ref(y);
  arena_matrix<promote_scalar_t<double, T>> y_val = value_of(y_ref);
  arena_matrix<plain_type_t<T>> arena_y = y_ref;

  const double r = y_val.norm();
  const double r_cubed = r * r * r;
  arena_matrix<plain_type_t<T>> res = y_val / r;

  reverse_pass_callback([arena_y, res, r, r_cubed, y_val]() mutable {
    const auto& adj = to_ref(res.adj());

    arena_y.adj()
        += adj / r - y_val * (y_val.array() * adj.array()).sum() / r_cubed;
  });

  return res;
}

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 * @param lp Log probability reference to increment.
 **/
template <int R, int C>
Eigen::Matrix<var, R, C> unit_vector_constrain(
    const Eigen::Matrix<var, R, C>& y, var& lp) {
  Eigen::Matrix<var, R, C> x = unit_vector_constrain(y);
  lp -= 0.5 * dot_self(y);
  return x;
}

}  // namespace math
}  // namespace stan
#endif
