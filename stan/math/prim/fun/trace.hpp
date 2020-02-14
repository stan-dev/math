#ifndef STAN_MATH_PRIM_FUN_TRACE_HPP
#define STAN_MATH_PRIM_FUN_TRACE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

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
template <typename T, typename = require_eigen_t<T>>
inline value_type_t<T> trace(const T& m) {
  return m.trace();
}

template <typename T, typename = require_not_eigen_t<T>, typename = void>
inline T trace(const T& m) {
  return m;
}

}  // namespace math
}  // namespace stan

#endif
