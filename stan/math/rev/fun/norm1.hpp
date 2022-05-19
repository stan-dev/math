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
 * @tparam EigVec type of the vector (must have one compile-time dimension equal to
 * 1)
 * @param[in] x Vector.
 * @return L1 norm of x.
 */
template <typename EigVec, require_eigen_vector_vt<is_var, EigVec>* = nullptr>
inline var norm1(const EigVec& x) {
  arena_t<EigVec> arena_v = x;
  var res = norm1(arena_v.val());
  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj().array() += res.adj() * sign(arena_v.val().array());
  });
  return res;
}

/**
 * Returns the L1 norm of a `var_value<Vector>`.
 *
 * @tparam VarVec `var_value<>` whose inner type has one compile-time row or column.
 * @param[in] x Vector.
 * @return L1 norm of x.
 */
//
template <typename VarVec, require_var_matrix_t<VarVec>* = nullptr>
inline var norm1(const VarVec& x) {
  return make_callback_vari(norm1(x.val()), [x](const auto& res) mutable {
    x.adj().array() += res.adj() * sign(x.val().array());
  });
}

}  // namespace math
}  // namespace stan
#endif
