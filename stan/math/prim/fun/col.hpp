#ifndef STAN_MATH_PRIM_FUN_COL_HPP
#define STAN_MATH_PRIM_FUN_COL_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the specified column of the specified matrix
 * using start-at-1 indexing.
 *
 * This is equivalent to calling <code>m.col(i - 1)</code> and
 * assigning the resulting template expression to a column vector.
 *
 * @tparam T type of elements in the matrix
 * @param m Matrix.
 * @param j Column index (count from 1).
 * @return Specified column of the matrix.
 * @throw std::out_of_range if j is out of range.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> col(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m, size_t j) {
  check_column_index("col", "j", m, j);
  return m.col(j - 1);
}

}  // namespace math
}  // namespace stan

#endif
