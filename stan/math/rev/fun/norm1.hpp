#ifndef STAN_MATH_REV_FUN_NORM1_HPP
#define STAN_MATH_REV_FUN_NORM1_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the L1 norm of a vector of var.
 *
 * @tparam T type of the vector (must have one compile-time dimension equal to
 * 1)
 * @param[in] v Vector.
 * @return L1 norm of v.
 */
template <typename T, require_eigen_vector_vt<is_var, T>* = nullptr>
inline var norm1(const T& v) {
  const auto& v_ref = to_ref(v);
  arena_t<T> arena_v(v_ref.size());
  double res_val = 0;
  for (size_t i = 0; i < arena_v.size(); ++i) {
    arena_v.coeffRef(i) = v_ref.coeffRef(i);
    res_val += abs(arena_v.coeffRef(i).val());
  }
  var res(res_val);

  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj().array() += res.adj() * sign(arena_v.val().array());
  });

  return res;
}

/**
 * Returns the L1 norm of a `var_value<Vector>`.
 *
 * @tparam A `var_value<>` whose inner type has one compile-time row or column.
 * @param[in] v Vector.
 * @return L1 norm of v.
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var norm1(const T& v) {
  var res = v.val().array().abs().sum();

  reverse_pass_callback([res, v]() mutable {
    v.adj().array() += res.adj() * sign(v.val().array());
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
