#ifndef STAN_MATH_REV_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
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
inline auto unit_vector_constrain(const T& y) {
  using ret_type = plain_type_t<T>;

  check_vector("unit_vector", "y", y);
  check_nonzero_size("unit_vector", "y", y);

  arena_t<T> arena_y = y;
  arena_t<T> y_val = arena_y.val();

  const double r = y_val.norm();
  arena_t<ret_type> res = y_val / r;

  reverse_pass_callback([arena_y, res, r, r_cubed, y_val]() mutable {
    arena_y.adj()
        += res.adj() / r
           - y_val * ((y_val.array() * res.adj().array()).sum() / (r * r * r));
  });

  return ret_type(res);
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
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto unit_vector_constrain(const T& y, var& lp) {
  const auto& y_ref = to_ref(y);
  auto x = unit_vector_constrain(y_ref);
  lp -= 0.5 * dot_self(y_ref);
  return x;
}

}  // namespace math
}  // namespace stan
#endif
