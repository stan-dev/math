#ifndef STAN_MATH_PRIM_FUN_COLS_HPP
#define STAN_MATH_PRIM_FUN_COLS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the number of columns in the specified
 * matrix, vector, or row vector.
 *
 * @tparam T type of the matrix
 * @param[in] m Input matrix, vector, or row vector.
 * @return Number of columns.
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline Eigen::Index cols(const T& m) {
  return m.cols();
}

}  // namespace math
}  // namespace stan

#endif
