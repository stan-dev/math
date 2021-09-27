#ifndef STAN_MATH_PRIM_FUN_ROWS_HPP
#define STAN_MATH_PRIM_FUN_ROWS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the number of rows in the specified
 * matrix, vector, or row vector.
 *
 * @tparam T type of the matrix
 * @param[in] m Input matrix, vector, or row vector.
 * @return Number of rows.
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline int rows(const T& m) {
  return m.rows();
}

}  // namespace math
}  // namespace stan

#endif
