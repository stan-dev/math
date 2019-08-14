#ifndef STAN_MATH_PRIM_MAT_FUN_TRACE_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRACE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * Returns the trace of the specified matrix.  The trace
 * is defined as the sum of the elements on the diagonal.
 * The matrix is not required to be square.  Returns 0 if
 * matrix is empty.
 *
 * @param[in] m Specified matrix.
 * @return Trace of the matrix.
 */
template <typename T, typename = enable_if_eigen<T>>
inline auto trace(const T& m) {
  return m.trace();
}

template <typename T, typename = enable_if_not_eigen<T>, typename = enable_if_arithmetic<T>>
inline auto trace(const T& m) {
  return m;
}
}  // namespace math
}  // namespace stan
#endif
