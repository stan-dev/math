#ifndef STAN_MATH_PRIM_MAT_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/log.hpp>
#include <stan/math/prim/mat/fun/sum.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log absolute determinant of the specified square matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m specified matrix
 * @return log absolute determinant of the matrix
 * @throw std::domain_error if matrix is not square and symmetric
 */
template <typename T, int R, int C>
inline T log_determinant_spd(const Eigen::Matrix<T, R, C>& m) {
  using std::log;
  check_symmetric("log_determinant_spd", "m", m);
  if (m.size() == 0)
    return 0;

  return sum(log(m.ldlt().vectorD().array()));
}

}  // namespace math
}  // namespace stan

#endif
