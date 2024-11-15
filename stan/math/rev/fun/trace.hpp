#ifndef STAN_MATH_REV_FUN_TRACE_HPP
#define STAN_MATH_REV_FUN_TRACE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/trace.hpp>

namespace stan {
namespace math {

/**
 * Returns the trace of the specified matrix.  The trace
 * is defined as the sum of the elements on the diagonal.
 * The matrix is not required to be square.  Returns 0 if
 * matrix is empty.
 *
 * @tparam T type of the elements in the matrix
 * @param[in] m Specified matrix.
 * @return Trace of the matrix.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto trace(const T& m) {
  arena_t<T> arena_m = m;

  return make_callback_var(arena_m.val_op().trace(),
                           [arena_m](const auto& vi) mutable {
                             arena_m.adj().diagonal().array() += vi.adj();
                           });
}

}  // namespace math
}  // namespace stan

#endif
