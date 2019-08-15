#ifndef STAN_MATH_PRIM_MAT_FUN_ROWS_DOT_SELF_HPP
#define STAN_MATH_PRIM_MAT_FUN_ROWS_DOT_SELF_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each row of a matrix with itself.
 * @param x Matrix.
 * @tparam T scalar type
 */
template <typename T, typename = enable_if_eigen<T>>
inline auto rows_dot_self(const T& x) {
  return x.rowwise().squaredNorm();
}

}  // namespace math
}  // namespace stan
#endif
