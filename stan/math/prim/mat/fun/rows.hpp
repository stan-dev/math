#ifndef STAN_MATH_PRIM_MAT_FUN_ROWS_HPP
#define STAN_MATH_PRIM_MAT_FUN_ROWS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the number of rows in the specified
 * matrix, vector, or row vector.
 *
 * @tparam T Type of matrix entries.
 * @tparam R Row type of matrix.
 * @tparam C Column type of matrix.
 * @param[in] m Input matrix, vector, or row vector.
 * @return Number of rows.
 */
template <typename T, typename = enable_if_eigen<T>>
inline auto rows(const T& m) {
  return m.rows();
}

}  // namespace math
}  // namespace stan
#endif
