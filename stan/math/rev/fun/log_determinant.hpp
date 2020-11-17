#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline var log_determinant(const T& m) {
  check_square("log_determinant", "m", m);

  if (m.size() == 0) {
    return var(0.0);
  }

  arena_t<T> arena_m = m;
  auto m_hh = arena_m.val().colPivHouseholderQr();
  auto arena_m_inv_transpose = to_arena(m_hh.inverse().transpose());
  var log_det = m_hh.logAbsDeterminant();

  reverse_pass_callback([arena_m, log_det, arena_m_inv_transpose]() mutable {
    arena_m.adj() += log_det.adj() * arena_m_inv_transpose;
  });
  return log_det;
}

}  // namespace math
}  // namespace stan
#endif
