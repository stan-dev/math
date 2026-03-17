#ifndef STAN_MATH_REV_FUN_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_DOT_SELF_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/square.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the dot product of a vector of var with itself.
 *
 * @tparam T type of the vector (must have one compile time dimension equal to
 * 1)
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename T, require_eigen_vector_vt<is_var, T>* = nullptr>
inline var dot_self(const T& v) {
  const auto& v_ref = to_ref(v);
  arena_t<T> arena_v(v_ref.size());
  double res_val = 0;
  for (size_t i = 0; i < arena_v.size(); ++i) {
    arena_v.coeffRef(i) = v_ref.coeffRef(i);
    res_val += arena_v.coeffRef(i).val() * arena_v.coeffRef(i).val();
  }
  var res(res_val);
  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj() += 2.0 * res.adj() * arena_v.val();
  });

  return res;
}

/**
 * Returns the dot product of a `var_value<Vector>` with itself.
 *
 * @tparam A `var_value<>` whose inner type has one compile time row or column.
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var dot_self(const T& v) {
  var res = v.val().dot(v.val());
  reverse_pass_callback(
      [res, v]() mutable { v.adj() += (2.0 * res.adj()) * v.val(); });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
