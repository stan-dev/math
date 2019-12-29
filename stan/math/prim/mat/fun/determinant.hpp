#ifndef STAN_MATH_PRIM_MAT_FUN_DETERMINANT_HPP
#define STAN_MATH_PRIM_MAT_FUN_DETERMINANT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the determinant of the specified square matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m Specified matrix.
 * @return Determinant of the matrix.
 * @throw std::domain_error if matrix is not square.
 */
template <typename T, int R, int C>
inline T determinant(const Eigen::Matrix<T, R, C>& m) {
  check_square("determinant", "m", m);
  if (m.size() == 0) {
    return 1;
  }

  return m.determinant();
}

}  // namespace math
}  // namespace stan

#endif
