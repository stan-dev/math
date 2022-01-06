#ifndef STAN_MATH_REV_FUN_NORM2_HPP
#define STAN_MATH_REV_FUN_NORM2_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the L2 norm of a vector of var.
 *
 * @tparam T type of the vector (must have one compile-time dimension equal to
 * 1)
 * @param[in] v Vector.
 * @return L2 norm of v.
 */
template <typename T, require_eigen_vector_vt<is_var, T>* = nullptr>
inline var norm2(const T& v) {
  arena_t<T> arena_v = v;
  var res = norm2(arena_v.val());
  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj().array() += res.adj() * (arena_v.val().array() / res.val());
  });
  return res;
}

/**
 * Returns the L2 norm of a `var_value<Vector>`.
 *
 * @tparam A `var_value<>` whose inner type has one compile-time row or column.
 * @param[in] v Vector.
 * @return L2 norm of v.
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var norm2(const T& v) {
  return make_callback_vari(norm2(v.val()), [v](const auto& res) mutable {
    v.adj().array() += res.adj() * (v.val().array() / res.val());
  });
}

}  // namespace math
}  // namespace stan
#endif
