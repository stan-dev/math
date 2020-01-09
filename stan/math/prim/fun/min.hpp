#ifndef STAN_MATH_PRIM_MAT_FUN_MIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_MIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the minimum coefficient in the specified
 * column vector.
 *
 * @param x specified vector
 * @return minimum coefficient value in the vector
 * @throws <code>std::invalid_argument</code> if the vector is size zero
 */
inline int min(const std::vector<int>& x) {
  check_nonzero_size("min", "int vector", x);
  Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 1>> m(&x[0], x.size());
  return m.minCoeff();
}

/**
 * Returns the minimum coefficient in the specified
 * column vector.
 *
 * @tparam T type of elements in the vector
 * @param x specified vector
 * @return minimum coefficient value in the vector, or infinity if the vector is
 * size zero
 */
template <typename T>
inline T min(const std::vector<T>& x) {
  if (x.size() == 0) {
    return INFTY;
  }
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> m(&x[0], x.size());
  return m.minCoeff();
}

/**
 * Returns the minimum coefficient in the specified
 * matrix, vector, or row vector.
 *
 * @tparam T type of elements in the matrix, vector or row vector
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param m specified matrix, vector, or row vector
 * @return minimum coefficient value in the vector, or infinity if the vector is
 * size zero
 */
template <typename T, int R, int C>
inline T min(const Eigen::Matrix<T, R, C>& m) {
  if (m.size() == 0) {
    return INFTY;
  }
  return m.minCoeff();
}

}  // namespace math
}  // namespace stan

#endif
