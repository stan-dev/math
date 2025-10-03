#ifndef STAN_MATH_REV_FUN_NORM2_HPP
#define STAN_MATH_REV_FUN_NORM2_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/norm2.hpp>

namespace stan {
namespace math {

/**
 * Returns the L2 norm of a vector of var.
 *
 * @tparam Container type of the vector (must have one compile-time dimension
 * equal to 1)
 * @param[in] x Vector.
 * @return L2 norm of x.
 */
template <typename Container,
          require_eigen_vector_vt<is_var, Container>* = nullptr>
inline var norm2(Container&& x) {
  arena_t<Container> arena_x = std::forward<Container>(x);
  var res = norm2(arena_x.val());
  reverse_pass_callback([res, arena_x]() mutable {
    arena_x.adj().array() += res.adj() * (arena_x.val().array() / res.val());
  });
  return res;
}

/**
 * Returns the L2 norm of a `var_value<Vector>`.
 *
 * @tparam Container `var_value<>` whose inner type has one compile-time row or
 * column.
 * @param[in] x Vector.
 * @return L2 norm of x.
 */
template <typename Container, require_var_matrix_t<Container>* = nullptr>
inline var norm2(const Container& x) {
  return make_callback_vari(norm2(x.val()), [x](const auto& res) mutable {
    x.adj().array() += res.adj() * (x.val().array() / res.val());
  });
}

}  // namespace math
}  // namespace stan
#endif
