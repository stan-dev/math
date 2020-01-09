#ifndef STAN_MATH_PRIM_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_PRIM_FUN_LOG_DETERMINANT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the log absolute determinant of the specified square matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m Specified matrix.
 * @return log absolute determinant of the matrix.
 * @throw std::domain_error if matrix is not square.
 */
template <typename T, int R, int C>
inline T log_determinant(const Eigen::Matrix<T, R, C>& m) {
  check_square("log_determinant", "m", m);
  if (m.size() == 0)
    return 0;

  return m.colPivHouseholderQr().logAbsDeterminant();
}

}  // namespace math
}  // namespace stan

#endif
