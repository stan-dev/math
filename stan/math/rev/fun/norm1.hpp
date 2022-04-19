#ifndef STAN_MATH_REV_FUN_NORM1_HPP
#define STAN_MATH_REV_FUN_NORM1_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sign.hpp>

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
  arena_t<T> arena_v = v;
  var res = norm1(arena_v.val());
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
//
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var norm1(const T& v) {
  return make_callback_vari(norm1(v.val()), [v](const auto& res) mutable {
    v.adj().array() += res.adj() * sign(v.val().array());
  });
}

}  // namespace math
}  // namespace stan
#endif
