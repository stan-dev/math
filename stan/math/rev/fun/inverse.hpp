#ifndef STAN_MATH_REV_FUN_INVERSE_HPP
#define STAN_MATH_REV_FUN_INVERSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Reverse mode specialization of calculating the inverse of the matrix.
 *
 * @param m specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not square.
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto inverse(const T& m) {
  using ret_type = plain_type_t<T>;
  
  check_square("inverse", "m", m);

  if (m.size() == 0) {
    return ret_type(0, 0);
  }

  arena_t<T> arena_m = m;
  auto arena_m_inv = to_arena(value_of(arena_m).inverse());
  arena_t<ret_type> res = arena_m_inv;

  reverse_pass_callback([res, arena_m_inv, arena_m]() mutable {
    arena_m.adj() -=
      arena_m_inv.transpose() * res.adj_op() * arena_m_inv.transpose();
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
