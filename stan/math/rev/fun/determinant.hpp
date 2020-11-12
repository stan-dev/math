#ifndef STAN_MATH_REV_FUN_DETERMINANT_HPP
#define STAN_MATH_REV_FUN_DETERMINANT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline var determinant(const T& m) {
  check_square("determinant", "m", m);

  if (m.size() == 0) {
    return var(1.0);
  }

  arena_t<T> arena_m = m;
  auto m_lu = arena_m.val().partialPivLu();
  auto arena_m_inv_transpose = to_arena(m_lu.inverse().transpose());
  var det = m_lu.determinant();

  reverse_pass_callback([arena_m, det, arena_m_inv_transpose]() mutable {
    arena_m.adj() += (det.adj() * det.val()) * arena_m_inv_transpose;
  });

  return det;
}

}  // namespace math
}  // namespace stan
#endif
